#!/usr/bin/env python3
"""
TAGE 预测器对比脚本
支持两种输入方式：
  1. CSV 文件：从 docs/example_and_reference_predictor_results.csv 过滤指定预测器
  2. 结果目录：读取 .out 文件

Usage:
  # 对比 CSV 中的 tage 和目录中的 my_bp
  python3 compare_predictors.py --base csv --base-predictor tage --compare dir --compare-dir results_my_bp/

  # 对比两个目录
  python3 compare_predictors.py --base dir --base-dir results/ --compare dir --compare-dir results_my_bp/

  # 对比 CSV 中的两个预测器
  python3 compare_predictors.py --base csv --base-predictor tage --compare csv --compare-predictor my_bp
"""

import sys
import os
import math
import csv
import argparse
from collections import defaultdict

# 预测器配置
PREDICTORS = {
    'tage': {
        'LOGLB': 6, 'NUMG': 8, 'LOGG': 11, 'LOGB': 12,
        'TAGW': 11, 'GHIST': 100, 'LOGP1': 14, 'GHIST1': 6
    },
    'my_bp': {
        'LOGLB': 6, 'NUMG': 8, 'LOGG': 12, 'LOGB': 12,
        'TAGW': 11, 'GHIST': 100, 'LOGP1': 14, 'GHIST1': 6
    }
}

MISPREDICTION_PENALTY = 8
DEFAULT_CSV = 'docs/example_and_reference_predictor_results.csv'

def calc_storage(params):
    """计算预测器存储开销"""
    LOGLB, NUMG, LOGG, LOGB, TAGW, GHIST, LOGP1, GHIST1 = (
        params['LOGLB'], params['NUMG'], params['LOGG'], params['LOGB'],
        params['TAGW'], params['GHIST'], params['LOGP1'], params['GHIST1']
    )

    LOGLINEINST = LOGLB - 2
    LINEINST = 1 << LOGLINEINST
    index1_bits = LOGP1 - LOGLINEINST
    bindex_bits = LOGB - LOGLINEINST
    HTAGBITS = TAGW - LOGLINEINST

    # 各组件存储计算 (bits)
    p1_total = (1 << index1_bits) * 1 * LINEINST
    gtag_total = (1 << LOGG) * TAGW * NUMG
    gpred_total = (1 << LOGG) * 1 * NUMG
    ghyst_total = (1 << LOGG) * 2 * NUMG
    ubit_total = (1 << LOGG) * 1 * NUMG
    bim_total = (1 << bindex_bits) * 1 * LINEINST

    return {
        'p1': p1_total,
        'gtag': gtag_total,
        'gpred': gpred_total,
        'ghyst': ghyst_total,
        'ubit': ubit_total,
        'bim': bim_total,
        'tage_subtotal': gtag_total + gpred_total + ghyst_total + ubit_total + bim_total,
        'total': p1_total + gtag_total + gpred_total + ghyst_total + ubit_total + bim_total,
        'params': {
            'LINEINST': LINEINST, 'index1_bits': index1_bits,
            'bindex_bits': bindex_bits, 'HTAGBITS': HTAGBITS
        }
    }

def parse_result_file(filepath):
    """解析 .out 文件"""
    try:
        with open(filepath, 'r') as f:
            line = f.readline().strip()
            if not line:
                return None
            parts = line.split(',')
            if len(parts) != 12:
                return None
            return {
                'name': parts[0],
                'instructions': int(parts[1]),
                'branches': int(parts[2]),
                'condbr': int(parts[3]),
                'npred': int(parts[4]),
                'extra': int(parts[5]),
                'diverge': int(parts[6]),
                'diverge_at_end': int(parts[7]),
                'misp': int(parts[8]),
                'p1_lat': float(parts[9]),
                'p2_lat': float(parts[10]),
                'epi': float(parts[11])
            }
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return None

def load_results_from_dir(directory):
    """从目录加载所有 .out 文件"""
    results = {}
    if not os.path.isdir(directory):
        print(f"Warning: Directory {directory} does not exist")
        return results
    for filename in os.listdir(directory):
        if not filename.endswith('.out'):
            continue
        filepath = os.path.join(directory, filename)
        data = parse_result_file(filepath)
        if data:
            results[data['name']] = data
    return results

def load_results_from_csv(csv_path, predictor_name):
    """从 CSV 文件加载指定预测器的结果"""
    results = {}
    if not os.path.isfile(csv_path):
        print(f"Error: CSV file {csv_path} does not exist")
        return results

    # CSV 列名映射（支持不同格式）
    col_map = {
        'trace': ['trace', 'trace_name', 'name', 'benchmark'],
        'instructions': ['instructions', 'instr'],
        'branches': ['branches', 'branch'],
        'condbr': ['conditional_branches', 'condbr'],
        'npred': ['predictions', 'npred'],
        'extra': ['extra_cycles', 'extra'],
        'diverge': ['p1_p2_disagreements', 'diverge'],
        'diverge_at_end': ['p1_p2_disagreements_at_end', 'diverge_at_end'],
        'misp': ['mispredictions', 'misp'],
        'p1_lat': ['p1_latency', 'p1_lat'],
        'p2_lat': ['p2_latency', 'p2_lat'],
        'epi': ['energy_per_instruction', 'epi', 'energy']
    }

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames

        # 构建列名映射
        def find_col(name):
            for h in headers:
                if h.lower() in col_map.get(name, []):
                    return h
            return name if name in headers else None

        predictor_col = 'predictor' if 'predictor' in headers else None

        for row in reader:
            # 过滤预测器
            if predictor_col and row.get(predictor_col, '') != predictor_name:
                continue

            trace_col = find_col('trace')
            if not trace_col or not row.get(trace_col):
                continue

            try:
                results[row[trace_col]] = {
                    'name': row[trace_col],
                    'instructions': int(row.get(find_col('instructions'), 0)),
                    'branches': int(row.get(find_col('branches'), 0)),
                    'condbr': int(row.get(find_col('condbr'), 0)),
                    'npred': int(row.get(find_col('npred'), 0)),
                    'extra': int(row.get(find_col('extra'), 0)),
                    'diverge': int(row.get(find_col('diverge'), 0)),
                    'diverge_at_end': int(row.get(find_col('diverge_at_end'), 0)),
                    'misp': int(row.get(find_col('misp'), 0)),
                    'p1_lat': float(row.get(find_col('p1_lat'), 0)),
                    'p2_lat': float(row.get(find_col('p2_lat'), 0)),
                    'epi': float(row.get(find_col('epi'), 0))
                }
            except (ValueError, KeyError) as e:
                print(f"Warning: Error parsing row for {row.get(trace_col, 'unknown')}: {e}")
                continue

    return results

def calc_metrics(data, p1_lat, p2_lat):
    """计算单个 benchmark 的指标"""
    if data['condbr'] == 0:
        accuracy = 1.0
    else:
        accuracy = 1.0 - (data['misp'] / data['condbr'])

    mpki = (data['misp'] / data['instructions']) * 1000 if data['instructions'] > 0 else 0
    misp_rate = data['misp'] / data['condbr'] if data['condbr'] > 0 else 0
    epi = data['epi']

    # 计算 IPC, CPI
    instructions = float(data['instructions'])
    pred_cycles = float(data['npred'])
    extra_cycles = float(data['extra'])
    divergences = float(data['diverge'])
    divergences_at_end = float(data['diverge_at_end'])
    mispredictions = float(data['misp'])

    if p2_lat <= p1_lat:
        cycles = pred_cycles * max(1, p2_lat)
    else:
        cycles = pred_cycles * max(1, p1_lat) + divergences * p2_lat - divergences_at_end * max(1, p1_lat)
    cycles += extra_cycles

    ipc = instructions / cycles if cycles > 0 else 0
    mpi = mispredictions / instructions if instructions > 0 else 0
    cpi = mpi * (MISPREDICTION_PENALTY + p2_lat)
    dpi = divergences / instructions if instructions > 0 else 0
    ppi = pred_cycles / instructions if instructions > 0 else 0

    return {
        'accuracy': accuracy,
        'mpki': mpki,
        'misp_rate': misp_rate,
        'epi': epi,
        'ipc': ipc,
        'cpi': cpi,
        'mpi': mpi,
        'dpi': dpi,
        'ppi': ppi
    }

def calc_vfs(ipc, cpi, epi):
    """计算 VFS 分数"""
    IPCcbp0 = 8
    CPIcbp0 = 0.0315
    EPIcbp0 = 1000
    ALPHA = 1.625
    BETA = 4 * ALPHA / (ALPHA - 1) ** 2
    GAMMA = 2 / (ALPHA - 1)
    cbp_energy_ratio = 0.05

    WPI0 = IPCcbp0 * CPIcbp0
    WPI = ipc * cpi
    speedup = (ipc / IPCcbp0) * (1 + WPI0) / (1 + WPI)
    LAMBDA = 1 / (1 + WPI0 / 2) - cbp_energy_ratio
    normalizedEPI = ((epi / EPIcbp0) * cbp_energy_ratio + LAMBDA * speedup ** GAMMA) * (1 + WPI / 2)

    if speedup * normalizedEPI <= 0:
        return 0

    vfs = speedup * ALPHA * (1 - 2 / (1 + math.sqrt(1 + BETA / (speedup * normalizedEPI))))
    return vfs

def print_storage_comparison(storage1, storage2, name1='tage', name2='my_bp'):
    """打印存储开销对比"""
    print("=" * 70)
    print("=== Storage Overhead Comparison ===")
    print("=" * 70)
    print(f"\n{name1} parameters: LOGLB={PREDICTORS[name1]['LOGLB']}, NUMG={PREDICTORS[name1]['NUMG']}, "
          f"LOGG={PREDICTORS[name1]['LOGG']}, LOGB={PREDICTORS[name1]['LOGB']}, "
          f"TAGW={PREDICTORS[name1]['TAGW']}, LOGP1={PREDICTORS[name1]['LOGP1']}")
    print(f"  LINEINST={storage1['params']['LINEINST']}, index1_bits={storage1['params']['index1_bits']}, "
          f"bindex_bits={storage1['params']['bindex_bits']}, HTAGBITS={storage1['params']['HTAGBITS']}")

    print(f"\n{name2} parameters: LOGLB={PREDICTORS[name2]['LOGLB']}, NUMG={PREDICTORS[name2]['NUMG']}, "
          f"LOGG={PREDICTORS[name2]['LOGG']}, LOGB={PREDICTORS[name2]['LOGB']}, "
          f"TAGW={PREDICTORS[name2]['TAGW']}, LOGP1={PREDICTORS[name2]['LOGP1']}")
    print(f"  LINEINST={storage2['params']['LINEINST']}, index1_bits={storage2['params']['index1_bits']}, "
          f"bindex_bits={storage2['params']['bindex_bits']}, HTAGBITS={storage2['params']['HTAGBITS']}")

    print("\n" + "-" * 70)
    print(f"{'Component':<18} {name1:>12} {name2:>12} {'Difference':>15}")
    print("-" * 70)

    comp_names = {
        'gtag': 'gtag', 'gpred': 'gpred', 'ghyst': 'ghyst',
        'ubit': 'ubit', 'bim': 'bim', 'p1': 'P1'
    }

    for key in ['gtag', 'gpred', 'ghyst', 'ubit', 'bim', 'p1']:
        v1 = storage1[key] / (1024 * 8)  # KB
        v2 = storage2[key] / (1024 * 8)
        diff_pct = ((v2 - v1) / v1 * 100) if v1 > 0 else 0
        diff_str = f"+{diff_pct:.1f}%" if diff_pct > 0 else f"{diff_pct:.1f}%"
        print(f"{comp_names[key]:<18} {v1:>10.2f} KB {v2:>10.2f} KB {diff_str:>15}")

    print("-" * 70)
    v1_total = storage1['total'] / (1024 * 8)
    v2_total = storage2['total'] / (1024 * 8)
    diff_total = ((v2_total - v1_total) / v1_total * 100) if v1_total > 0 else 0
    diff_str = f"+{diff_total:.1f}%" if diff_total > 0 else f"{diff_total:.1f}%"
    print(f"{'Total':<18} {v1_total:>10.2f} KB {v2_total:>10.2f} KB {diff_str:>15}")
    print("=" * 70)

def print_per_benchmark_comparison(common_benchmarks, metrics1, metrics2,
                                   data1, data2, name1='tage', name2='my_bp'):
    """打印每个 benchmark 的对比"""
    print("\n" + "=" * 130)
    print("=== Per-Benchmark Comparison ===")
    print("=" * 130)

    header = f"{'Benchmark':<30} {name1+' Acc':>11} {name2+' Acc':>11} {'Δ Acc':>9} "
    header += f"{name1+' MPKI':>11} {name2+' MPKI':>11} {'Δ MPKI':>9} "
    header += f"{name1+' EPI':>9} {name2+' EPI':>9} {'Norm':>7} "
    header += f"{name1+' IPC':>9} {name2+' IPC':>9} {'Δ IPC':>9} "
    header += f"{name1+' VFS':>9} {name2+' VFS':>9} {'Δ VFS':>9}"
    print(header)
    print("-" * 130)

    rows = []
    for bench in sorted(common_benchmarks):
        m1 = metrics1[bench]
        m2 = metrics2[bench]

        d_acc = (m2['accuracy'] - m1['accuracy']) * 100
        d_mpki = ((m2['mpki'] - m1['mpki']) / m1['mpki'] * 100) if m1['mpki'] > 0 else 0
        norm_epi = m2['epi'] / m1['epi'] if m1['epi'] > 0 else 0
        d_ipc = ((m2['ipc'] - m1['ipc']) / m1['ipc'] * 100) if m1['ipc'] > 0 else 0

        v1 = calc_vfs(m1['ipc'], m1['cpi'], m1['epi'])
        v2 = calc_vfs(m2['ipc'], m2['cpi'], m2['epi'])
        d_vfs = ((v2 - v1) / v1 * 100) if v1 > 0 else 0

        row = {
            'benchmark': bench,
            'acc1': m1['accuracy'], 'acc2': m2['accuracy'], 'd_acc': d_acc,
            'mpki1': m1['mpki'], 'mpki2': m2['mpki'], 'd_mpki': d_mpki,
            'epi1': m1['epi'], 'epi2': m2['epi'], 'norm_epi': norm_epi,
            'ipc1': m1['ipc'], 'ipc2': m2['ipc'], 'd_ipc': d_ipc,
            'vfs1': v1, 'vfs2': v2, 'd_vfs': d_vfs,
            'instr': data1[bench]['instructions'],
            'condbr': data1[bench]['condbr'],
            'misp1': data1[bench]['misp'], 'misp2': data2[bench]['misp']
        }
        rows.append(row)

    # 按 accuracy 差异排序，显示改进最大的在前
    rows.sort(key=lambda x: x['d_acc'], reverse=True)

    for row in rows:
        line = f"{row['benchmark']:<30} {row['acc1']*100:>9.2f}% {row['acc2']*100:>9.2f}% "
        line += f"{'+' if row['d_acc']>0 else ''}{row['d_acc']:.3f}% "
        line += f"{row['mpki1']:>9.3f} {row['mpki2']:>9.3f} "
        line += f"{'+' if row['d_mpki']>0 else ''}{row['d_mpki']:.2f}% "
        line += f"{row['epi1']:>7.1f} {row['epi2']:>7.1f} {row['norm_epi']:>7.4f} "
        line += f"{row['ipc1']:>7.4f} {row['ipc2']:>7.4f} "
        line += f"{'+' if row['d_ipc']>0 else ''}{row['d_ipc']:.2f}% "
        line += f"{row['vfs1']:>7.4f} {row['vfs2']:>7.4f} "
        line += f"{'+' if row['d_vfs']>0 else ''}{row['d_vfs']:.2f}%"
        print(line)

    print("=" * 130)
    return rows

def print_overall_statistics(common_benchmarks, metrics1, metrics2, name1='tage', name2='my_bp'):
    """打印总体统计"""
    print("\n" + "=" * 70)
    print("=== Overall Statistics ===")
    print("=" * 70)

    # 收集所有指标
    acc1_list, acc2_list = [], []
    mpki1_list, mpki2_list = [], []
    epi1_list, epi2_list = [], []
    ipc1_list, ipc2_list = [], []
    cpi1_list, cpi2_list = [], []
    vfs1_list, vfs2_list = [], []

    for bench in common_benchmarks:
        m1 = metrics1[bench]
        m2 = metrics2[bench]
        acc1_list.append(m1['accuracy'])
        acc2_list.append(m2['accuracy'])
        mpki1_list.append(m1['mpki'])
        mpki2_list.append(m2['mpki'])
        epi1_list.append(m1['epi'])
        epi2_list.append(m2['epi'])
        ipc1_list.append(m1['ipc'])
        ipc2_list.append(m2['ipc'])
        cpi1_list.append(m1['cpi'])
        cpi2_list.append(m2['cpi'])

        v1 = calc_vfs(m1['ipc'], m1['cpi'], m1['epi'])
        v2 = calc_vfs(m2['ipc'], m2['cpi'], m2['epi'])
        vfs1_list.append(v1)
        vfs2_list.append(v2)

    # 几何平均准确率
    geo_mean_acc1 = math.exp(sum(math.log(a) for a in acc1_list) / len(acc1_list))
    geo_mean_acc2 = math.exp(sum(math.log(a) for a in acc2_list) / len(acc2_list))

    # 算术平均
    avg_mpki1 = sum(mpki1_list) / len(mpki1_list)
    avg_mpki2 = sum(mpki2_list) / len(mpki2_list)
    avg_epi1 = sum(epi1_list) / len(epi1_list)
    avg_epi2 = sum(epi2_list) / len(epi2_list)

    # 调和平均 IPC
    harmonic_ipc1 = len(ipc1_list) / sum(1/i for i in ipc1_list)
    harmonic_ipc2 = len(ipc2_list) / sum(1/i for i in ipc2_list)

    avg_cpi1 = sum(cpi1_list) / len(cpi1_list)
    avg_cpi2 = sum(cpi2_list) / len(cpi2_list)

    avg_vfs1 = sum(vfs1_list) / len(vfs1_list)
    avg_vfs2 = sum(vfs2_list) / len(vfs2_list)

    print(f"\n{'Metric':<25} {name1:>15} {name2:>15} {'Relative Improvement':>20}")
    print("-" * 70)

    d_acc = (geo_mean_acc2 - geo_mean_acc1) * 100
    print(f"Geometric Mean Accuracy {geo_mean_acc1*100:>13.4f}% {geo_mean_acc2*100:>13.4f}% "
          f"{'+' if d_acc>0 else ''}{d_acc:.4f}%")

    d_mpki = ((avg_mpki2 - avg_mpki1) / avg_mpki1 * 100) if avg_mpki1 > 0 else 0
    print(f"Average MPKI              {avg_mpki1:>13.4f} {avg_mpki2:>13.4f} "
          f"{'+' if d_mpki>0 else ''}{d_mpki:.2f}% (lower is better)")

    d_epi = ((avg_epi2 - avg_epi1) / avg_epi1 * 100) if avg_epi1 > 0 else 0
    print(f"Average EPI (fJ)          {avg_epi1:>13.2f} {avg_epi2:>13.2f} "
          f"{'+' if d_epi>0 else ''}{d_epi:.2f}%")

    d_ipc = ((harmonic_ipc2 - harmonic_ipc1) / harmonic_ipc1 * 100) if harmonic_ipc1 > 0 else 0
    print(f"Harmonic Mean IPC         {harmonic_ipc1:>13.4f} {harmonic_ipc2:>13.4f} "
          f"{'+' if d_ipc>0 else ''}{d_ipc:.2f}%")

    d_cpi = ((avg_cpi2 - avg_cpi1) / avg_cpi1 * 100) if avg_cpi1 > 0 else 0
    print(f"Average CPI               {avg_cpi1:>13.6f} {avg_cpi2:>13.6f} "
          f"{'+' if d_cpi>0 else ''}{d_cpi:.2f}% (lower is better)")

    d_vfs = ((avg_vfs2 - avg_vfs1) / avg_vfs1 * 100) if avg_vfs1 > 0 else 0
    print(f"Average VFS Score         {avg_vfs1:>13.4f} {avg_vfs2:>13.4f} "
          f"{'+' if d_vfs>0 else ''}{d_vfs:.2f}%")

    print("=" * 70)

    return {
        'geo_mean_acc1': geo_mean_acc1, 'geo_mean_acc2': geo_mean_acc2,
        'avg_mpki1': avg_mpki1, 'avg_mpki2': avg_mpki2,
        'avg_epi1': avg_epi1, 'avg_epi2': avg_epi2,
        'harmonic_ipc1': harmonic_ipc1, 'harmonic_ipc2': harmonic_ipc2,
        'avg_cpi1': avg_cpi1, 'avg_cpi2': avg_cpi2,
        'avg_vfs1': avg_vfs1, 'avg_vfs2': avg_vfs2
    }

def export_csv(common_benchmarks, metrics1, metrics2, data1, data2,
               output_file='comparison_report.csv', name1='tage', name2='my_bp'):
    """导出 CSV 格式报告"""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)

        # 写入表头
        header = ['Benchmark',
                  f'{name1}_accuracy', f'{name2}_accuracy', 'accuracy_diff',
                  f'{name1}_mpki', f'{name2}_mpki', 'mpki_diff_pct',
                  f'{name1}_epi', f'{name2}_epi', 'epi_ratio',
                  f'{name1}_ipc', f'{name2}_ipc', 'ipc_diff_pct',
                  f'{name1}_cpi', f'{name2}_cpi', 'cpi_diff_pct',
                  f'{name1}_vfs', f'{name2}_vfs', 'vfs_diff_pct',
                  'instructions', 'condbr',
                  f'{name1}_misp', f'{name2}_misp', 'misp_diff']
        writer.writerow(header)

        for bench in sorted(common_benchmarks):
            m1 = metrics1[bench]
            m2 = metrics2[bench]

            v1 = calc_vfs(m1['ipc'], m1['cpi'], m1['epi'])
            v2 = calc_vfs(m2['ipc'], m2['cpi'], m2['epi'])

            row = [
                bench,
                m1['accuracy'], m2['accuracy'], (m2['accuracy'] - m1['accuracy']),
                m1['mpki'], m2['mpki'], ((m2['mpki'] - m1['mpki']) / m1['mpki'] * 100) if m1['mpki'] > 0 else 0,
                m1['epi'], m2['epi'], (m2['epi'] / m1['epi']) if m1['epi'] > 0 else 0,
                m1['ipc'], m2['ipc'], ((m2['ipc'] - m1['ipc']) / m1['ipc'] * 100) if m1['ipc'] > 0 else 0,
                m1['cpi'], m2['cpi'], ((m2['cpi'] - m1['cpi']) / m1['cpi'] * 100) if m1['cpi'] > 0 else 0,
                v1, v2, ((v2 - v1) / v1 * 100) if v1 > 0 else 0,
                data1[bench]['instructions'],
                data1[bench]['condbr'],
                data1[bench]['misp'],
                data2[bench]['misp'],
                data2[bench]['misp'] - data1[bench]['misp']
            ]
            writer.writerow(row)

    print(f"\nCSV report exported to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Compare two branch predictors')

    # Base (reference) predictor
    parser.add_argument('--base', choices=['csv', 'dir'], required=True,
                        help='Source type for base predictor')
    parser.add_argument('--base-csv', default=DEFAULT_CSV,
                        help='Path to CSV file (default: docs/example_and_reference_predictor_results.csv)')
    parser.add_argument('--base-predictor', default='tage',
                        help='Predictor name to filter from CSV (default: tage)')
    parser.add_argument('--base-dir', default='results/',
                        help='Directory with .out files for base predictor')

    # Compare predictor
    parser.add_argument('--compare', choices=['csv', 'dir'], required=True,
                        help='Source type for compare predictor')
    parser.add_argument('--compare-csv', default=DEFAULT_CSV,
                        help='Path to CSV file for compare predictor')
    parser.add_argument('--compare-predictor', default='my_bp',
                        help='Predictor name to filter from CSV (default: my_bp)')
    parser.add_argument('--compare-dir', default='results_my_bp/',
                        help='Directory with .out files for compare predictor')

    # Output
    parser.add_argument('--output', default='comparison_report.csv',
                        help='Output CSV file path (default: comparison_report.csv)')
    parser.add_argument('--name1', default='tage',
                        help='Display name for base predictor (default: tage)')
    parser.add_argument('--name2', default='my_bp',
                        help='Display name for compare predictor (default: my_bp)')

    args = parser.parse_args()

    # Load base results
    print(f"Loading base predictor '{args.base_predictor}' results...")
    if args.base == 'csv':
        results1 = load_results_from_csv(args.base_csv, args.base_predictor)
        print(f"  Loaded from CSV: {args.base_csv}")
    else:
        results1 = load_results_from_dir(args.base_dir)
        print(f"  Loaded from directory: {args.base_dir}")
    print(f"  Total: {len(results1)} benchmarks")

    # Load compare results
    print(f"\nLoading compare predictor '{args.compare_predictor}' results...")
    if args.compare == 'csv':
        results2 = load_results_from_csv(args.compare_csv, args.compare_predictor)
        print(f"  Loaded from CSV: {args.compare_csv}")
    else:
        results2 = load_results_from_dir(args.compare_dir)
        print(f"  Loaded from directory: {args.compare_dir}")
    print(f"  Total: {len(results2)} benchmarks")

    # Find common benchmarks
    common = set(results1.keys()) & set(results2.keys())
    print(f"\nCommon benchmarks: {len(common)}")

    if len(common) == 0:
        print("ERROR: No common benchmarks found!")
        print(f"  Base has: {list(results1.keys())[:5]}...")
        print(f"  Compare has: {list(results2.keys())[:5]}...")
        sys.exit(1)

    name1 = args.name1
    name2 = args.name2

    # 如果两个预测器都在 PREDICTORS 配置中，打印存储对比
    if name1 in PREDICTORS and name2 in PREDICTORS:
        storage1 = calc_storage(PREDICTORS[name1])
        storage2 = calc_storage(PREDICTORS[name2])
        print_storage_comparison(storage1, storage2, name1, name2)
    else:
        print("\nNote: Storage comparison skipped (predictor configs not found)")

    # 计算 P1/P2 延迟
    p1_lat = 0
    p2_lat = 0
    for bench in common:
        p1_lat = max(p1_lat, math.ceil(results1[bench]['p1_lat']))
        p1_lat = max(p1_lat, math.ceil(results2[bench]['p1_lat']))
        p2_lat = max(p2_lat, math.ceil(results1[bench]['p2_lat']))
        p2_lat = max(p2_lat, math.ceil(results2[bench]['p2_lat']))

    print(f"\nUsing latencies: P1={p1_lat}, P2={p2_lat}")

    # 计算每个 benchmark 的指标
    metrics1 = {b: calc_metrics(results1[b], p1_lat, p2_lat) for b in common}
    metrics2 = {b: calc_metrics(results2[b], p1_lat, p2_lat) for b in common}

    # 打印每 benchmark 对比
    print_per_benchmark_comparison(common, metrics1, metrics2, results1, results2, name1, name2)

    # 打印总体统计
    overall = print_overall_statistics(common, metrics1, metrics2, name1, name2)

    # 导出 CSV
    export_csv(common, metrics1, metrics2, results1, results2, args.output, name1, name2)

    print("\nDone!")

if __name__ == '__main__':
    main()

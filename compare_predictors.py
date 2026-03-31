#!/usr/bin/env python3
"""
TAGE 预测器对比脚本
支持两种输入方式：
  1. CSV 文件：从 docs/example_and_reference_predictor_results.csv 过滤指定预测器
  2. 结果目录：读取 .out 文件

特性：
  - 对比 base 和 compare 两个预测器的性能
  - 自动加载 reference MPKI 作为参考基准
  - 支持同名预测器修改前后对比（如修改 LOGG 参数）
  - 支持手动指定存储配置参数

使用场景:

# 场景 1: 对比 tage 和 my_bp (不同预测器)
python3 compare_predictors.py --base csv --base-predictor tage --compare dir --compare-dir results_my_bp/

# 场景 2: 对比同名预测器修改前后
python3 compare_predictors.py \\
  --base dir --base-dir results_tage_v1/ \\
  --compare dir --compare-dir results_tage_v2/ \\
  --base-label "LOGG=11" --compare-label "LOGG=12" \\
  --storage1 "LOGLB=6,NUMG=8,LOGG=11,LOGB=12,TAGW=11" \\
  --storage2 "LOGLB=6,NUMG=8,LOGG=12,LOGB=12,TAGW=11"

# 场景 3: 对比 tage 和 reference
python3 compare_predictors.py --base csv --base-predictor tage --compare csv --compare-predictor reference

# 场景 4: 对比两个目录的结果
python3 compare_predictors.py --base dir --base-dir results/ --compare dir --compare-dir results_my_bp/
"""

import sys
import os
import math
import csv
import argparse
import re
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

DEFAULT_CSV = 'docs/example_and_reference_predictor_results.csv'
P2_TO_EXEC_STAGES = 9

VFS_REF_IPC = 8.0
VFS_REF_CPI = 0.0315
VFS_REF_EPI = 1000.0
VFS_ALPHA = 1.625
VFS_BETA = 4 * VFS_ALPHA / (VFS_ALPHA - 1) ** 2
VFS_GAMMA = 2 / (VFS_ALPHA - 1)
VFS_CBP_ENERGY_RATIO = 0.05

def load_ref_mpki_from_csv(csv_path):
    """从 CSV 加载 reference 预测器的 MPKI"""
    ref_mpki = {}
    if not os.path.isfile(csv_path):
        return ref_mpki

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get('predictor', '') != 'reference':
                continue
            trace = row.get('trace', '')
            if not trace:
                continue
            try:
                misp = int(row.get('mispredictions', 0))
                instr = int(row.get('instructions', 0))
                ref_mpki[trace] = (misp / instr * 1000) if instr > 0 else 0
            except (ValueError, KeyError):
                continue
    return ref_mpki

def parse_storage_config(config_str):
    """解析存储配置字符串，如 'LOGLB=6,NUMG=8,LOGG=11'"""
    if not config_str:
        return None
    params = {}
    for item in config_str.split(','):
        if '=' in item:
            key, val = item.split('=', 1)
            params[key.strip()] = int(val.strip())
    return params

def infer_trace_group(trace_name):
    """将 trace 名称映射到组名（含 5-prefix 特例）"""
    # 用户要求：5 开头 trace 归为同一组
    if trace_name.startswith('5'):
        return 'spec5xx'

    group = trace_name

    # step1: 去除末尾 "_digits"
    if '_' in group:
        prefix, suffix = group.rsplit('_', 1)
        if suffix.isdigit() and prefix:
            group = prefix

    # step2: 去除末尾 ".digits"
    group = re.sub(r'\.\d+$', '', group)
    # step3: 去除末尾 "-digits"
    group = re.sub(r'-\d+$', '', group)

    return group if group else trace_name

def build_trace_groups(benchmarks):
    """根据 trace 名称构建分组"""
    groups = defaultdict(list)
    for bench in benchmarks:
        groups[infer_trace_group(bench)].append(bench)
    for traces in groups.values():
        traces.sort()
    return dict(sorted(groups.items(), key=lambda kv: kv[0]))

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
        'config': dict(params),
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

def compute_predictor_level_latencies(results, benchmarks):
    """计算目录级 P1/P2 最大延迟（不取整，保留浮点）"""
    p1_lat = 0.0
    p2_lat = 0.0
    for bench in benchmarks:
        p1_lat = max(p1_lat, float(results[bench]['p1_lat']))
        p2_lat = max(p2_lat, float(results[bench]['p2_lat']))
    return p1_lat, p2_lat

def compute_trace_level_latencies(data):
    """按单条 trace 的输出延迟（浮点原值）"""
    return float(data['p1_lat']), float(data['p2_lat'])

def calc_prediction_cycles(data, p1_lat, p2_lat):
    """按 predictor_metrics.py 当前公式计算预测周期数"""
    pred_cycles = float(data['npred'])
    extra_cycles = float(data['extra'])
    divergences = float(data['diverge'])
    divergences_at_end = float(data['diverge_at_end'])

    if p2_lat <= p1_lat:
        cycles = pred_cycles * max(1, p2_lat)
    else:
        cycles = (
            pred_cycles * max(1, p1_lat)
            + divergences * p2_lat
            - divergences_at_end * max(1, p1_lat)
        )

    return cycles + extra_cycles

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
    mispredictions = float(data['misp'])

    cycles = calc_prediction_cycles(data, p1_lat, p2_lat)

    ipc = instructions / cycles if cycles > 0 else 0
    mpi = mispredictions / instructions if instructions > 0 else 0
    cpi = mpi * (P2_TO_EXEC_STAGES + p2_lat - max(1, min(p1_lat, p2_lat)))
    divergences = float(data['diverge'])
    pred_cycles = float(data['npred'])
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

def calc_metrics_per_trace_latency(data):
    """按该 trace 自身的 P1/P2 延迟计算指标"""
    trace_p1_lat, trace_p2_lat = compute_trace_level_latencies(data)
    return calc_metrics(data, trace_p1_lat, trace_p2_lat)

def aggregate_results(results, benchmarks, name):
    """将多条 trace 的原始计数器聚合为一条记录"""
    agg = {
        'name': name,
        'instructions': 0,
        'branches': 0,
        'condbr': 0,
        'npred': 0,
        'extra': 0,
        'diverge': 0,
        'diverge_at_end': 0,
        'misp': 0,
        'p1_lat': 0.0,
        'p2_lat': 0.0,
        'epi': 0.0,
    }
    total_instructions = 0
    weighted_epi = 0.0

    for bench in benchmarks:
        d = results[bench]
        agg['instructions'] += int(d['instructions'])
        agg['branches'] += int(d['branches'])
        agg['condbr'] += int(d['condbr'])
        agg['npred'] += int(d['npred'])
        agg['extra'] += int(d['extra'])
        agg['diverge'] += int(d['diverge'])
        agg['diverge_at_end'] += int(d['diverge_at_end'])
        agg['misp'] += int(d['misp'])
        agg['p1_lat'] = max(agg['p1_lat'], float(d['p1_lat']))
        agg['p2_lat'] = max(agg['p2_lat'], float(d['p2_lat']))

        instr = int(d['instructions'])
        total_instructions += instr
        weighted_epi += float(d['epi']) * instr

    agg['epi'] = (weighted_epi / total_instructions) if total_instructions > 0 else 0.0
    return agg

def aggregate_ref_mpki(ref_mpki, benchmarks, instr_results):
    """用指令数加权聚合 reference MPKI"""
    ref_misp = 0.0
    ref_instr = 0.0
    for bench in benchmarks:
        if bench not in ref_mpki:
            continue
        instr = float(instr_results[bench]['instructions'])
        ref_misp += (float(ref_mpki[bench]) / 1000.0) * instr
        ref_instr += instr
    if ref_instr <= 0:
        return None
    return (ref_misp / ref_instr) * 1000.0

def summarize_predictor_metrics(results, benchmarks):
    """按 per-trace 延迟口径汇总 common benchmarks 上的指标"""
    p1_lat, p2_lat = compute_predictor_level_latencies(results, benchmarks)

    count = 0
    sum_inverse_ipc = 0.0
    sum_cpi = 0.0
    sum_epi = 0.0
    sum_mpi = 0.0
    sum_dpi = 0.0
    sum_ppi = 0.0

    for bench in benchmarks:
        metrics = calc_metrics_per_trace_latency(results[bench])
        if metrics['ipc'] <= 0:
            continue
        count += 1
        sum_inverse_ipc += 1.0 / metrics['ipc']
        sum_cpi += metrics['cpi']
        sum_epi += metrics['epi']
        sum_mpi += metrics['mpi']
        sum_dpi += metrics['dpi']
        sum_ppi += metrics['ppi']

    if count == 0:
        return {
            'ipc': 0.0,
            'cpi': 0.0,
            'epi': 0.0,
            'mpi': 0.0,
            'dpi': 0.0,
            'ppi': 0.0,
            'mpki': 0.0,
            'p1_latency': p1_lat,
            'p2_latency': p2_lat,
            'vfs': 0.0,
        }

    ipc = count / sum_inverse_ipc
    cpi = sum_cpi / count
    epi = sum_epi / count
    mpi = sum_mpi / count
    dpi = sum_dpi / count
    ppi = sum_ppi / count

    return {
        'ipc': ipc,
        'cpi': cpi,
        'epi': epi,
        'mpi': mpi,
        'dpi': dpi,
        'ppi': ppi,
        'mpki': mpi * 1000.0,
        'p1_latency': p1_lat,
        'p2_latency': p2_lat,
        'vfs': calc_vfs(ipc, cpi, epi),
    }

def calc_vfs(ipc, cpi, epi):
    """计算 VFS 分数"""
    WPI0 = VFS_REF_IPC * VFS_REF_CPI
    WPI = ipc * cpi
    speedup = (ipc / VFS_REF_IPC) * (1 + WPI0) / (1 + WPI)
    LAMBDA = 1 / (1 + WPI0 / 2) - VFS_CBP_ENERGY_RATIO
    normalizedEPI = (
        ((epi / VFS_REF_EPI) * VFS_CBP_ENERGY_RATIO + LAMBDA * speedup ** VFS_GAMMA)
        * (1 + WPI / 2)
    )

    if speedup * normalizedEPI <= 0:
        return 0

    vfs = speedup * VFS_ALPHA * (1 - 2 / (1 + math.sqrt(1 + VFS_BETA / (speedup * normalizedEPI))))
    return vfs

def format_storage_config(storage):
    """格式化打印存储配置"""
    config = storage.get('config', {})
    ordered_keys = ['LOGLB', 'NUMG', 'LOGG', 'LOGB', 'TAGW', 'GHIST', 'LOGP1', 'GHIST1']
    parts = [f"{key}={config[key]}" for key in ordered_keys if key in config]
    return ', '.join(parts) if parts else 'N/A'

def print_storage_comparison(storage1, storage2, name1='tage', name2='my_bp'):
    """打印存储开销对比"""
    print("=" * 70)
    print("=== Storage Overhead Comparison ===")
    print("=" * 70)
    print(f"\n{name1} parameters: {format_storage_config(storage1)}")
    print(f"  LINEINST={storage1['params']['LINEINST']}, index1_bits={storage1['params']['index1_bits']}, "
          f"bindex_bits={storage1['params']['bindex_bits']}, HTAGBITS={storage1['params']['HTAGBITS']}")

    print(f"\n{name2} parameters: {format_storage_config(storage2)}")
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
                                   data1, data2, ref_mpki,
                                   name1='tage', name2='my_bp'):
    """打印每个 benchmark 的对比"""
    print("\n" + "=" * 155)
    print("=== Per-Benchmark Comparison ===")
    print("=" * 155)

    header = f"{'Benchmark':<28} {'Ref MPKI':>10} {name1+' MPKI':>11} {name2+' MPKI':>11} "
    header += f"{'vs Ref1':>9} {'vs Ref2':>9} "
    header += f"{name1+' Acc':>10} {name2+' Acc':>10} {'Δ Acc':>9} "
    header += f"{name1+' EPI':>8} {name2+' EPI':>8} {'Norm':>7} "
    header += f"{name1+' VFS':>8} {name2+' VFS':>8} {'Δ VFS':>9}"
    print(header)
    print("-" * 155)

    rows = []
    for bench in sorted(common_benchmarks):
        m1 = metrics1[bench]
        m2 = metrics2[bench]
        ref = ref_mpki.get(bench, None)

        d_acc = (m2['accuracy'] - m1['accuracy']) * 100
        norm_epi = m2['epi'] / m1['epi'] if m1['epi'] > 0 else 0

        # vs reference MPKI
        d_mpki1_ref = ((m1['mpki'] - ref) / ref * 100) if ref and ref > 0 else 0
        d_mpki2_ref = ((m2['mpki'] - ref) / ref * 100) if ref and ref > 0 else 0

        v1 = calc_vfs(m1['ipc'], m1['cpi'], m1['epi'])
        v2 = calc_vfs(m2['ipc'], m2['cpi'], m2['epi'])
        d_vfs = ((v2 - v1) / v1 * 100) if v1 > 0 else 0

        row = {
            'benchmark': bench,
            'ref_mpki': ref,
            'mpki1': m1['mpki'], 'mpki2': m2['mpki'],
            'd_mpki1_ref': d_mpki1_ref, 'd_mpki2_ref': d_mpki2_ref,
            'acc1': m1['accuracy'], 'acc2': m2['accuracy'], 'd_acc': d_acc,
            'epi1': m1['epi'], 'epi2': m2['epi'], 'norm_epi': norm_epi,
            'vfs1': v1, 'vfs2': v2, 'd_vfs': d_vfs,
            'instr': data1[bench]['instructions'],
            'condbr': data1[bench]['condbr'],
            'misp1': data1[bench]['misp'], 'misp2': data2[bench]['misp']
        }
        rows.append(row)

    # 按 accuracy 差异排序，显示改进最大的在前
    rows.sort(key=lambda x: x['d_acc'], reverse=True)

    for row in rows:
        ref_str = f"{row['ref_mpki']:.3f}" if row['ref_mpki'] is not None else "N/A"
        line = f"{row['benchmark']:<28} {ref_str:>10} {row['mpki1']:>9.3f} {row['mpki2']:>9.3f} "
        line += f"{'+' if row['d_mpki1_ref']>0 else ''}{row['d_mpki1_ref']:>7.1f}% "
        line += f"{'+' if row['d_mpki2_ref']>0 else ''}{row['d_mpki2_ref']:>7.1f}% "
        line += f"{row['acc1']*100:>8.2f}% {row['acc2']*100:>8.2f}% "
        line += f"{'+' if row['d_acc']>0 else ''}{row['d_acc']:>7.3f}% "
        line += f"{row['epi1']:>6.1f} {row['epi2']:>6.1f} {row['norm_epi']:>7.4f} "
        line += f"{row['vfs1']:>6.4f} {row['vfs2']:>6.4f} "
        line += f"{'+' if row['d_vfs']>0 else ''}{row['d_vfs']:>7.2f}%"
        print(line)

    print("=" * 155)
    return rows

def print_per_group_comparison(group_members, results1, results2, ref_mpki,
                               name1='tage', name2='my_bp'):
    """打印每个 trace 组的聚合对比"""
    print("\n" + "=" * 170)
    print("=== Per-Group Comparison (Aggregated by Trace Prefix) ===")
    print("=" * 170)

    header = f"{'Group':<28} {'#Traces':>8} {'Ref MPKI':>10} {name1+' MPKI':>11} {name2+' MPKI':>11} "
    header += f"{'Δ MPKI':>9} {name1+' Acc':>10} {name2+' Acc':>10} {'Δ Acc':>9} "
    header += f"{name1+' EPI':>9} {name2+' EPI':>9} {name1+' VFS':>9} {name2+' VFS':>9} {'Δ VFS':>9}"
    print(header)
    print("-" * 170)

    rows = []
    for group, traces in group_members.items():
        gdata1 = aggregate_results(results1, traces, group)
        gdata2 = aggregate_results(results2, traces, group)

        gp1_1, gp2_1 = compute_predictor_level_latencies(results1, traces)
        gp1_2, gp2_2 = compute_predictor_level_latencies(results2, traces)
        m1 = calc_metrics(gdata1, gp1_1, gp2_1)
        m2 = calc_metrics(gdata2, gp1_2, gp2_2)
        v1 = calc_vfs(m1['ipc'], m1['cpi'], m1['epi'])
        v2 = calc_vfs(m2['ipc'], m2['cpi'], m2['epi'])

        d_mpki = ((m2['mpki'] - m1['mpki']) / m1['mpki'] * 100) if m1['mpki'] > 0 else 0.0
        d_acc = (m2['accuracy'] - m1['accuracy']) * 100.0
        d_vfs = ((v2 - v1) / v1 * 100) if v1 > 0 else 0.0
        g_ref_mpki = aggregate_ref_mpki(ref_mpki, traces, results1)

        row = {
            'group': group,
            'num_traces': len(traces),
            'ref_mpki': g_ref_mpki,
            'mpki1': m1['mpki'],
            'mpki2': m2['mpki'],
            'd_mpki': d_mpki,
            'acc1': m1['accuracy'],
            'acc2': m2['accuracy'],
            'd_acc': d_acc,
            'epi1': m1['epi'],
            'epi2': m2['epi'],
            'vfs1': v1,
            'vfs2': v2,
            'd_vfs': d_vfs,
            'instructions': gdata1['instructions'],
            'condbr': gdata1['condbr'],
            'misp1': gdata1['misp'],
            'misp2': gdata2['misp'],
            'members': ';'.join(traces),
        }
        rows.append(row)

    rows.sort(key=lambda x: x['d_acc'], reverse=True)
    for row in rows:
        ref_str = f"{row['ref_mpki']:.3f}" if row['ref_mpki'] is not None else "N/A"
        line = f"{row['group']:<28} {row['num_traces']:>8} {ref_str:>10} "
        line += f"{row['mpki1']:>9.3f} {row['mpki2']:>9.3f} "
        line += f"{'+' if row['d_mpki'] > 0 else ''}{row['d_mpki']:>7.2f}% "
        line += f"{row['acc1']*100:>8.2f}% {row['acc2']*100:>8.2f}% "
        line += f"{'+' if row['d_acc'] > 0 else ''}{row['d_acc']:>7.3f}% "
        line += f"{row['epi1']:>7.1f} {row['epi2']:>7.1f} "
        line += f"{row['vfs1']:>7.4f} {row['vfs2']:>7.4f} "
        line += f"{'+' if row['d_vfs'] > 0 else ''}{row['d_vfs']:>7.2f}%"
        print(line)

    print("=" * 170)
    return rows

def print_overall_statistics(common_benchmarks, metrics1, metrics2, ref_mpki,
                             summary1, summary2,
                             name1='tage', name2='my_bp'):
    """打印总体统计"""
    print("\n" + "=" * 70)
    print("=== Overall Statistics ===")
    print("=" * 70)

    # 收集所有指标
    acc1_list, acc2_list = [], []
    mpki_ref_list = []

    for bench in common_benchmarks:
        m1 = metrics1[bench]
        m2 = metrics2[bench]
        ref = ref_mpki.get(bench, None)

        acc1_list.append(m1['accuracy'])
        acc2_list.append(m2['accuracy'])
        if ref is not None:
            mpki_ref_list.append(ref)

    # 几何平均准确率
    geo_mean_acc1 = math.exp(sum(math.log(a) for a in acc1_list) / len(acc1_list))
    geo_mean_acc2 = math.exp(sum(math.log(a) for a in acc2_list) / len(acc2_list))

    avg_mpki1 = summary1['mpki']
    avg_mpki2 = summary2['mpki']
    avg_mpki_ref = sum(mpki_ref_list) / len(mpki_ref_list) if mpki_ref_list else None

    avg_epi1 = summary1['epi']
    avg_epi2 = summary2['epi']
    harmonic_ipc1 = summary1['ipc']
    harmonic_ipc2 = summary2['ipc']
    avg_cpi1 = summary1['cpi']
    avg_cpi2 = summary2['cpi']
    avg_vfs1 = summary1['vfs']
    avg_vfs2 = summary2['vfs']
    avg_dpi1 = summary1['dpi']
    avg_dpi2 = summary2['dpi']
    avg_ppi1 = summary1['ppi']
    avg_ppi2 = summary2['ppi']

    print(f"\n{'Metric':<25} {name1:>15} {name2:>15} {'vs Ref':>10} {'vs Ref':>10}")
    print("-" * 70)

    d_acc = (geo_mean_acc2 - geo_mean_acc1) * 100
    print(f"Geometric Mean Accuracy {geo_mean_acc1*100:>13.4f}% {geo_mean_acc2*100:>13.4f}% "
          f"{'+' if d_acc>0 else ''}{d_acc:.4f}%")

    if avg_mpki_ref is not None:
        d_mpki1_ref = ((avg_mpki1 - avg_mpki_ref) / avg_mpki_ref * 100) if avg_mpki_ref > 0 else 0
        d_mpki2_ref = ((avg_mpki2 - avg_mpki_ref) / avg_mpki_ref * 100) if avg_mpki_ref > 0 else 0
        print(f"Average MPKI              {avg_mpki1:>13.4f} {avg_mpki2:>13.4f} "
              f"{'+' if d_mpki1_ref>0 else ''}{d_mpki1_ref:.2f}% "
              f"{'+' if d_mpki2_ref>0 else ''}{d_mpki2_ref:.2f}%")
    else:
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

    d_dpi = ((avg_dpi2 - avg_dpi1) / avg_dpi1 * 100) if avg_dpi1 > 0 else 0
    print(f"Average DPI               {avg_dpi1:>13.6f} {avg_dpi2:>13.6f} "
          f"{'+' if d_dpi>0 else ''}{d_dpi:.2f}% (lower is better)")

    d_ppi = ((avg_ppi2 - avg_ppi1) / avg_ppi1 * 100) if avg_ppi1 > 0 else 0
    print(f"Average PPI               {avg_ppi1:>13.6f} {avg_ppi2:>13.6f} "
          f"{'+' if d_ppi>0 else ''}{d_ppi:.2f}% (lower is better)")

    print(f"P1/P2 Latency (cycles)   {summary1['p1_latency']:>6}/{summary1['p2_latency']:<6} "
          f"{summary2['p1_latency']:>6}/{summary2['p2_latency']:<6}")

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
        'avg_dpi1': avg_dpi1, 'avg_dpi2': avg_dpi2,
        'avg_ppi1': avg_ppi1, 'avg_ppi2': avg_ppi2,
        'avg_vfs1': avg_vfs1, 'avg_vfs2': avg_vfs2
    }

def export_csv(common_benchmarks, metrics1, metrics2, data1, data2, ref_mpki,
               output_file='comparison_report.csv', name1='tage', name2='my_bp'):
    """导出 CSV 格式报告"""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)

        # 写入表头
        header = ['Benchmark',
                  'ref_mpki',
                  f'{name1}_accuracy', f'{name2}_accuracy', 'accuracy_diff',
                  f'{name1}_mpki', f'{name2}_mpki', 'mpki_diff_pct',
                  f'{name1}_mpki_vs_ref', f'{name2}_mpki_vs_ref',
                  f'{name1}_epi', f'{name2}_epi', 'epi_ratio',
                  f'{name1}_vfs', f'{name2}_vfs', 'vfs_diff_pct',
                  'instructions', 'condbr',
                  f'{name1}_misp', f'{name2}_misp', 'misp_diff']
        writer.writerow(header)

        for bench in sorted(common_benchmarks):
            m1 = metrics1[bench]
            m2 = metrics2[bench]
            ref = ref_mpki.get(bench, None)

            v1 = calc_vfs(m1['ipc'], m1['cpi'], m1['epi'])
            v2 = calc_vfs(m2['ipc'], m2['cpi'], m2['epi'])

            # vs ref MPKI
            mpki1_vs_ref = ((m1['mpki'] - ref) / ref * 100) if ref and ref > 0 else None
            mpki2_vs_ref = ((m2['mpki'] - ref) / ref * 100) if ref and ref > 0 else None

            row = [
                bench,
                ref if ref is not None else '',
                m1['accuracy'], m2['accuracy'], (m2['accuracy'] - m1['accuracy']),
                m1['mpki'], m2['mpki'], ((m2['mpki'] - m1['mpki']) / m1['mpki'] * 100) if m1['mpki'] > 0 else 0,
                mpki1_vs_ref if mpki1_vs_ref is not None else '',
                mpki2_vs_ref if mpki2_vs_ref is not None else '',
                m1['epi'], m2['epi'], (m2['epi'] / m1['epi']) if m1['epi'] > 0 else 0,
                v1, v2, ((v2 - v1) / v1 * 100) if v1 > 0 else 0,
                data1[bench]['instructions'],
                data1[bench]['condbr'],
                data1[bench]['misp'],
                data2[bench]['misp'],
                data2[bench]['misp'] - data1[bench]['misp']
            ]
            writer.writerow(row)

    print(f"\nCSV report exported to: {output_file}")

def derive_group_output_path(output_file):
    base, ext = os.path.splitext(output_file)
    if not ext:
        ext = '.csv'
    return f"{base}_groups{ext}"

def export_group_csv(group_rows, output_file='comparison_report_groups.csv',
                     name1='tage', name2='my_bp'):
    """导出分组聚合报告"""
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Group', 'num_traces', 'ref_mpki',
            f'{name1}_accuracy', f'{name2}_accuracy', 'accuracy_diff',
            f'{name1}_mpki', f'{name2}_mpki', 'mpki_diff_pct',
            f'{name1}_epi', f'{name2}_epi', 'epi_ratio',
            f'{name1}_vfs', f'{name2}_vfs', 'vfs_diff_pct',
            'instructions', 'condbr', f'{name1}_misp', f'{name2}_misp', 'misp_diff',
            'members'
        ])

        for row in sorted(group_rows, key=lambda x: x['group']):
            writer.writerow([
                row['group'],
                row['num_traces'],
                row['ref_mpki'] if row['ref_mpki'] is not None else '',
                row['acc1'], row['acc2'], row['d_acc'],
                row['mpki1'], row['mpki2'], row['d_mpki'],
                row['epi1'], row['epi2'], (row['epi2'] / row['epi1']) if row['epi1'] > 0 else 0.0,
                row['vfs1'], row['vfs2'], row['d_vfs'],
                row['instructions'],
                row['condbr'],
                row['misp1'],
                row['misp2'],
                row['misp2'] - row['misp1'],
                row['members']
            ])

    print(f"Grouped CSV report exported to: {output_file}")

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

    # 存储对比配置（支持同名预测器不同参数）
    parser.add_argument('--storage1', default=None,
                        help='Storage config for predictor 1 (format: LOGLB=6,NUMG=8,LOGG=11,...)')
    parser.add_argument('--storage2', default=None,
                        help='Storage config for predictor 2 (format: LOGLB=6,NUMG=8,LOGG=12,...)')

    # Output
    parser.add_argument('--output', default='comparison_report.csv',
                        help='Output CSV file path (default: comparison_report.csv)')
    parser.add_argument('--name1', default=None,
                        help='Display name for base predictor (default: auto from predictor param)')
    parser.add_argument('--name2', default=None,
                        help='Display name for compare predictor (default: auto from predictor param)')
    parser.add_argument('--base-label', default='base',
                        help='Label for base in output (e.g., "before", "v1", default: "base")')
    parser.add_argument('--compare-label', default='modified',
                        help='Label for compare in output (e.g., "after", "v2", default: "modified")')
    parser.add_argument('--reference-csv', default=DEFAULT_CSV,
                        help='Path to CSV file used to load reference MPKI (default: docs/example_and_reference_predictor_results.csv)')

    args = parser.parse_args()

    # 确定显示名称
    # 如果用户没有指定 name1/name2，根据场景自动选择
    if args.name1 is None:
        name1 = args.base_label if args.base_predictor == args.compare_predictor else args.base_predictor
    else:
        name1 = args.name1

    if args.name2 is None:
        name2 = args.compare_label if args.base_predictor == args.compare_predictor else args.compare_predictor
    else:
        name2 = args.name2

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

    group_members = build_trace_groups(common)
    print(f"Trace groups: {len(group_members)}")

    # 存储对比：优先使用用户指定的配置，否则从 PREDICTORS 获取
    storage1 = None
    storage2 = None

    if args.storage1:
        storage1 = calc_storage(parse_storage_config(args.storage1))
    elif args.base_predictor in PREDICTORS:
        storage1 = calc_storage(PREDICTORS[args.base_predictor])

    if args.storage2:
        storage2 = calc_storage(parse_storage_config(args.storage2))
    elif args.compare_predictor in PREDICTORS:
        storage2 = calc_storage(PREDICTORS[args.compare_predictor])

    if storage1 and storage2:
        print_storage_comparison(storage1, storage2, name1, name2)
    else:
        print("\nNote: Storage comparison skipped (predictor configs not found or use --storage1/--storage2)")

    # 加载 reference MPKI
    print(f"\nLoading reference MPKI from {args.reference_csv}...")
    ref_mpki = load_ref_mpki_from_csv(args.reference_csv)
    print(f"  Loaded {len(ref_mpki)} benchmarks")

    summary1 = summarize_predictor_metrics(results1, common)
    summary2 = summarize_predictor_metrics(results2, common)

    print(f"\n{name1} predictor_metrics latencies: P1={summary1['p1_latency']}, P2={summary1['p2_latency']}")
    print(f"{name2} predictor_metrics latencies: P1={summary2['p1_latency']}, P2={summary2['p2_latency']}")

    # 计算每个 benchmark 的指标
    metrics1 = {b: calc_metrics_per_trace_latency(results1[b]) for b in common}
    metrics2 = {b: calc_metrics_per_trace_latency(results2[b]) for b in common}

    # 打印每 benchmark 对比（包含 reference MPKI）
    print_per_benchmark_comparison(common, metrics1, metrics2, results1, results2, ref_mpki, name1, name2)

    # 打印每组聚合对比
    group_rows = print_per_group_comparison(group_members, results1, results2, ref_mpki, name1, name2)

    # 打印总体统计
    overall = print_overall_statistics(common, metrics1, metrics2, ref_mpki, summary1, summary2, name1, name2)

    # 导出 CSV（包含 reference MPKI）
    export_csv(common, metrics1, metrics2, results1, results2, ref_mpki, args.output, name1, name2)
    export_group_csv(group_rows, derive_group_output_path(args.output), name1, name2)

    print("\nDone!")

if __name__ == '__main__':
    main()

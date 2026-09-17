import pandas as pd
import os

def filter_continuous_boundaries(input_path, output_path, boundary_only_path, 
                                 insulation_threshold=1.4, 
                                 direction_threshold=0.3):
    df = pd.read_csv(input_path)
    
    # 确保列存在
    required_columns = ['Boundary', 'insulation_score', 'direction_score']
    for col in required_columns:
        if col not in df.columns:
            raise KeyError(f"Missing column: {col}")
    
    # 1. 基础硬性过滤（只保留生物学上合格的边界）
    condition = (df['insulation_score'] > insulation_threshold) & (df['direction_score'].abs() < direction_threshold)
    df_filtered = df[condition].sort_values(by='Boundary').reset_index(drop=True)
    
    if df_filtered.empty:
        print(f"无符合条件的TAD边界: {input_path}")
        return None
    
    print(f"初步过滤后剩余 {len(df_filtered)} 个候选边界点")
    
    # 2. 连续性分组（找到相邻编号差值为1的簇）
    groups = []
    current_group = [df_filtered.iloc[0]]
    for i in range(1, len(df_filtered)):
        if df_filtered.iloc[i]['Boundary'] == df_filtered.iloc[i-1]['Boundary'] + 1:
            current_group.append(df_filtered.iloc[i])
        else:
            groups.append(pd.DataFrame(current_group))
            current_group = [df_filtered.iloc[i]]
    groups.append(pd.DataFrame(current_group))
    
    print(f"共分为 {len(groups)} 个连续簇")
    
    # 3. 每个簇内：基于锐度（Sharpness）加权排序，选出最优边界
    best_boundaries = []
    for group in groups:
        if len(group) == 1:
            # 单个点直接保留
            best = group.iloc[0]
        else:
            # 组内有多个连续点，计算每个点的局部锐度（相对于组内最小分数）
            group_copy = group.copy()
            min_score = group_copy['insulation_score'].min()
            max_score = group_copy['insulation_score'].max()
            
            # 锐度计算：当前点相比组内最低分的提升幅度（值越大越尖锐）
            # 如果组内分数都一样，锐度设为 0
            if max_score - min_score > 1e-6:
                group_copy['sharpness'] = (group_copy['insulation_score'] - min_score) / (max_score - min_score)
            else:
                group_copy['sharpness'] = 0.0
            
            # 综合得分 = 绝缘分数 * (1 + 锐度系数)
            # 这样既保证分数高的优先，又确保尖锐的峰能排在平缓坡之前
            group_copy['combined_score'] = group_copy['insulation_score'] * (1 + group_copy['sharpness'])
            
            # 按综合得分降序排列，取第一名
            group_sorted = group_copy.sort_values(by='combined_score', ascending=False)
            best = group_sorted.iloc[0]
        
        best_boundaries.append(best)
    
    # 4. 构建最终结果
    result = pd.DataFrame(best_boundaries)[['Boundary', 'insulation_score', 'direction_score']]
    result = result.sort_values(by='Boundary').reset_index(drop=True)
    
    # 保存完整结果（含所有列）
    result.to_csv(output_path, index=False)
    
    # 输出边界ID列表（仅用于第三部分）
    with open(boundary_only_path, 'w') as f:
        for boundary in result['Boundary']:
            f.write(f"{int(boundary)}\n")
    
    print(f"筛选完成: {input_path}，最终保留 {len(result)} 个高精度边界")
    return result




if __name__ == "__main__":
    input_parent_dir = "...path.../results"         # ✅
    output_root = "...path.../results/filtered_results"   # ✅
    boundary_root = "...path.../results/boundary_only"   # ✅

    resolutions = [25, 50, 100]
    res_subdirs = [f"{res}k" for res in resolutions]  
    

    chromosomes = [20, 21, 22]
    
    for res, res_subdir in zip(resolutions, res_subdirs):
        print(f"\n===== Processing resolution: {res_subdir} =====")
        
        input_dir = os.path.join(input_parent_dir, res_subdir)
        if not os.path.exists(input_dir):
            print(f"Warning: Input directory for resolution {res_subdir} does not exist → {input_dir}, skipping this resolution")
            continue

        output_dir = os.path.join(output_root, res_subdir)
        boundary_dir = os.path.join(boundary_root, res_subdir)
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(boundary_dir, exist_ok=True)
        print(f"Output directory: {output_dir}, Boundary directory: {boundary_dir}")
        
        for chr_num in chromosomes:
            print(f"\n--- chr{chr_num} ---")
            
          
            input_path = os.path.join(input_dir, f"chr{chr_num}_clustering_results.txt")
          
            output_path = os.path.join(output_dir, f"chr{chr_num}_best_continuous_boundaries.txt")
            boundary_only_path = os.path.join(boundary_dir, f"chr{chr_num}_best_boundaries_only.txt")
            
          
            if not os.path.exists(input_path):
                print(f"Warning: Input file for chromosome chr{chr_num} does not exist → {input_path}, skipping")
                continue
            
          
            try:
                print("input_path:",input_path)
                print("output_path:",output_path)
                print("boundary_only_path:",boundary_only_path)
                filter_continuous_boundaries(input_path, output_path, boundary_only_path)
                print(f"Processing completed: {input_path} → output to {output_path} and {boundary_only_path}")
            except Exception as e:
                print(f"Error processing chromosome chr{chr_num}: {e}, skipping")
                continue
import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

def get_neighbors(node_1based, incidence_matrix, direction=None):
    """
    Parameters:
        node_1based: Current node (1-based)
        direction: 'upstream' (smaller node number)、'downstream' (larger node number)
    """
    neighbors = set()
   
    node_idx = node_1based - 1  
    for j in range(incidence_matrix.shape[1]):
        if incidence_matrix[node_idx, j] > 0:  
            edge_nodes_0based = np.where(incidence_matrix[:, j] > 0)[0]
            edge_nodes_1based = [n + 1 for n in edge_nodes_0based]  
            if direction == 'upstream':
                edge_nodes_1based = [n for n in edge_nodes_1based if n < node_1based]
            elif direction == 'downstream':
                edge_nodes_1based = [n for n in edge_nodes_1based if n > node_1based]
            
            neighbors.update(edge_nodes_1based)
    
    neighbors.discard(node_1based)
    return sorted(neighbors)

def calculate_shared_edges_weight(node_list_1based, incidence_matrix):
    """Calculate the sum of weights of shared hyperedges for all node pairs in the 1-based node list"""
    shared_edges_weight = 0
   
    for i in range(len(node_list_1based)):
        for j in range(i + 1, len(node_list_1based)):
            n1 = node_list_1based[i]
            n2 = node_list_1based[j]
            n1_idx = n1 - 1
            n2_idx = n2 - 1
            
            common_edges = np.where((incidence_matrix[n1_idx] > 0) & (incidence_matrix[n2_idx] > 0))[0]
            for edge in common_edges:
                shared_edges_weight += min(incidence_matrix[n1_idx, edge], incidence_matrix[n2_idx, edge])
    
    return shared_edges_weight

def compute_insulation_and_direction(boundary_node_1based, incidence_matrix, max_distance=40):
    """
    基于超图通量的绝缘分数与方向性计算
    
    Parameters:
        boundary_node_1based: 候选边界节点（1-based）
        incidence_matrix: 超图关联矩阵 (nodes x hyperedges)
        max_distance: 上下游邻居的最大距离限制
    
    Returns:
        dict: {'insulation_score': float, 'direction_score': float}
    """
    # 获取上下游邻居
    upstream_neighbors = get_neighbors(boundary_node_1based, incidence_matrix, direction='upstream')
    downstream_neighbors = get_neighbors(boundary_node_1based, incidence_matrix, direction='downstream')
    
    upstream_neighbors = [n for n in upstream_neighbors if abs(n - boundary_node_1based) <= max_distance]
    downstream_neighbors = [n for n in downstream_neighbors if abs(n - boundary_node_1based) <= max_distance]
    
    # 若某一侧无节点，则无法形成有效边界，直接返回最低分
    if not upstream_neighbors or not downstream_neighbors:
        return {
            'insulation_score': -999.0, 
            'direction_score': 0.0,
            'upstream_weight': 0.0,
            'downstream_weight': 0.0,
            'cross_weight': 0.0 
        }
    
    # 转换为0-based索引用于矩阵切片
    up_idx = [n - 1 for n in upstream_neighbors]
    down_idx = [n - 1 for n in downstream_neighbors]
    
    # 提取子矩阵
    sub_U = incidence_matrix[up_idx, :]  # shape: (len(up), E)
    sub_D = incidence_matrix[down_idx, :]  # shape: (len(down), E)
    
    # 按列求和，得到每条超边落在左右两侧的总权重
    sum_U = np.sum(sub_U, axis=0)  # 长度 = 超边总数 E
    sum_D = np.sum(sub_D, axis=0)
    
    # 核心分类：判断每条超边属于“左侧内部”、“右侧内部”还是“跨边界”
    mask_U_only = (sum_U > 0) & (sum_D == 0)
    mask_D_only = (sum_U == 0) & (sum_D > 0)
    mask_cross = (sum_U > 0) & (sum_D > 0)
    
    # 计算总权重
    weight_inside_U = np.sum(sum_U[mask_U_only])
    weight_inside_D = np.sum(sum_D[mask_D_only])
    weight_cross = np.sum(sum_U[mask_cross] + sum_D[mask_cross])  # 两侧权重相加，保留完整Pore-C信号
    
    # 对数绝缘分数（加1平滑，防止log(0)）
    insulation_score = np.log((weight_inside_U + weight_inside_D + 1.0) / (weight_cross + 1.0))
    
    # 方向性评分（用于过滤掉非TAD边界的噪音点）
    total_signal = weight_inside_U + weight_inside_D + weight_cross
    if total_signal > 0:
        direction_score = (weight_inside_D - weight_inside_U) / total_signal
    else:
        direction_score = 0.0
    
    return {
        'insulation_score': insulation_score,
        'direction_score': direction_score,
        'upstream_weight': weight_inside_U,
        'downstream_weight': weight_inside_D,
        'cross_weight': weight_cross
    }



def load_all_lines_as_matrix(file2):
    lines = []
    with open(file2, 'r', encoding='utf-8') as file:
        for line in file:
            lines.append([int(value) for value in line.strip().split()])
    matrix = np.array(lines)
    return matrix


def read_boundary_nodes(file_path):
   
    with open(file_path, 'r') as f:
        lines = f.readlines()
    boundaries = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        boundaries.append(int(line))  
    return boundaries

def save_results_to_file(results, output_path):
    """
    Format: Boundary, Upstream clustering coefficient, Cross clustering coefficient, Downstream clustering coefficient

    """
    with open(output_path, 'w') as f:
 
        f.write("Boundary,insulation_score,direction_score,cross_weight\n")

        for res in results:
            line = f"{res['boundary_node']},{res['insulation_score']},{res['direction_score']},{res['cross_weight']}\n"
            f.write(line)
    print(f"Results have been saved to the file: {output_path}")


def load_all_lines_as_matrix(file2):
    lines = []
    with open(file2, 'r', encoding='utf-8') as file:
        for line in file:
            lines.append([int(value) for value in line.strip().split()])
    return np.array(lines)


def read_boundary_nodes(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    boundaries = []
    for line in lines:
        line = line.strip()
        if line:
            boundaries.append(int(line))  
    return boundaries


def process_single_task(args):
   
    res, chr_num, base_matrix_path, base_boundary_prefix, output_dir, res_to_maxdistance = args
    try:
        pid = os.getpid()
        max_distance = res_to_maxdistance[res]
        res_tag = f"{res}k"
        print(f"[Process {pid}] Starting processing: resolution {res_tag}, chromosome chr{chr_num} maximum distance: {max_distance}")


        output_subdir = os.path.join(output_dir, res_tag)
        os.makedirs(output_subdir, exist_ok=True)
        
        boundary_res_dir = os.path.join(base_boundary_prefix, f"no_results_{res}")
        matrix_path = os.path.join(base_matrix_path, f"chr{chr_num}_{res}_VE_matrix.txt")
        boundary_path = os.path.join(boundary_res_dir, f"chr{chr_num}_boundary.txt")
        output_path = os.path.join(output_subdir, f"chr{chr_num}_clustering_results.txt")
     
        if not os.path.exists(boundary_res_dir):
            raise FileNotFoundError(f"Boundary directory does not exist: {boundary_res_dir}")
        if not os.path.exists(matrix_path):
            raise FileNotFoundError(f"Matrix file does not exist: {matrix_path}")
        if not os.path.exists(boundary_path):
            raise FileNotFoundError(f"Boundary file does not exist: {boundary_path}")
        
      
        incidence_matrix = load_all_lines_as_matrix(matrix_path)
        total_nodes = incidence_matrix.shape[0]
      
        boundary_nodes = read_boundary_nodes(boundary_path)
      

        valid_boundaries = [n for n in boundary_nodes if 1 <= n <= total_nodes]
        invalid_nodes = [n for n in boundary_nodes if not (1 <= n <= total_nodes)]
        if invalid_nodes:
            print(f"[Process {pid}] Warning: skipping invalid nodes: {invalid_nodes}")
        if not valid_boundaries:
            print(f"[Process {pid}] No valid boundary nodes, skipping")
            return (res, chr_num, "No valid boundary nodes")
        
        
        all_results = []
        
        for node in valid_boundaries:
            res1 = compute_insulation_and_direction(node, incidence_matrix, max_distance)
            result_dict = {
                'boundary_node': node,
                'insulation_score': res1['insulation_score'],   # 核心指标（越高越好）
                'direction_score': res1['direction_score'],     # 辅助筛选
                'cross_weight': res1['cross_weight']            # 保留用于调试
            }
            all_results.append(result_dict)
    
        
        print(f"Completed calculation of clustering coefficients for chromosome {chr_num}, with a total of {len(all_results)} boundary nodes")
        
  
        save_results_to_file(all_results, output_path)
    
    except Exception as e:
        print(f"[Process {os.getpid()}] Error processing resolution {res}k chromosome {chr_num}: {e}")
        return (res, chr_num, f"Failed: {str(e)}")


def main():
   
    base_matrix_path = "...path.../GM12878/hypergraph"     # ✅### 完整的超图加权关联矩阵

    base_boundary_prefix = "...path..."      # ✅  保存文件根目录 
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)

    res_to_maxdistance = {25: 80, 50: 40, 100: 20}
 
    chromosomes = [20,21,22]  #
    resolutions = [25,50,100] #,100
    max_outer_workers = 4  
    
 
    tasks = [
        (res, chr_num, base_matrix_path, base_boundary_prefix, output_dir, res_to_maxdistance)
        for res in resolutions
        for chr_num in chromosomes
    ]
    print(f"Total number of tasks: {len(tasks)}, running {max_outer_workers} tasks simultaneously")
 
    with ProcessPoolExecutor(max_workers=max_outer_workers) as outer_executor:
        futures = {outer_executor.submit(process_single_task, task): task for task in tasks}
  
        results = []
        for future in as_completed(futures):
            task = futures[future]
            res, chr_num = task[0], task[1]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                results.append((res, chr_num, f"Task submission failed: {str(e)}"))

    print("\n===== All tasks have been processed! =====")
    # for res, chr_num, status in results:
    #     print(f"Resolution {res}k, chromosome {chr_num}: {status}")


if __name__ == "__main__":
    main()
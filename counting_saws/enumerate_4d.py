import os
import json

def load_progress(file_path):
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            data = json.load(f)
            max_n = data['max_n']
            counts = {int(k): v for k, v in data['counts'].items()}
    else:
        max_n = 0
        counts = {}
    return max_n, counts

def save_progress(file_path, counts):
    max_n = max(counts.keys())
    with open(file_path, 'w') as f:
        json.dump({'max_n': max_n, 'counts': counts}, f)

def is_valid_move(position, visited_positions):
    return position not in visited_positions

def enumerate_walks(n, current_position, visited_positions, counts, current_length):
    if n == current_length:
        if current_length in counts:
            counts[current_length] += 1
        else:
            counts[current_length] = 1
        return
    
    directions = [
        (1, 0, 0, 0), (-1, 0, 0, 0),  # x-axis
        (0, 1, 0, 0), (0, -1, 0, 0),  # y-axis
        (0, 0, 1, 0), (0, 0, -1, 0),  # z-axis
        (0, 0, 0, 1), (0, 0, 0, -1)   # w-axis
    ]
    
    for direction in directions:
        next_position = (
            current_position[0] + direction[0],
            current_position[1] + direction[1],
            current_position[2] + direction[2],
            current_position[3] + direction[3]
        )
        if is_valid_move(next_position, visited_positions):
            visited_positions.add(next_position)
            enumerate_walks(n + 1, next_position, visited_positions, counts, current_length)
            visited_positions.remove(next_position)

def compute_saws(N_max):
    file_path = 'SAW_progress_4d.json'
    max_n, counts = load_progress(file_path)
    
    start_n = max_n + 1
    
    for n in range(start_n, N_max + 1):
        visited_positions = {(0, 0, 0, 0)}
        enumerate_walks(1, (0, 0, 0, 0), visited_positions, counts, n)
        save_progress(file_path, counts)
    
    print(f"Computed counts for N_max={N_max}: {counts}")
    return counts

if __name__ == "__main__":
    N_max = 6  # Change this value to set the desired maximum length of SAWs
    counts = compute_saws(N_max)

import os
import glob
import csv
import sys
import traceback

# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src'))

from run_learning_is_v2 import main as run_learning

def get_models():
    base_dir = "models/phase_1"
    models = []
    for root, dirs, files in os.walk(base_dir):
        for f in files:
            if f.endswith(".json"):
                models.append(os.path.join(root, f))
    return sorted(models)

def main():
    models = get_models()
    output_file = "benchmark_results_learning.csv"
    
    print(f"Found {len(models)} models.")
    
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['model', 'method', 'probability', 'std_error', 'ci_lower', 'ci_upper', 'total_time', 'status']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for model_path in models:
            model_name = os.path.basename(model_path)
            print("="*60)
            print(f"Processing Model: {model_name}")
            print(f"Path: {model_path}")
            
            method_name = "Learning-Based IS v2"
            print(f"  Running {method_name}...")
            
            try:
                res = run_learning(model_path)
                
                row = {
                    'model': model_name,
                    'method': method_name,
                    'probability': res.get('probability'),
                    'std_error': res.get('std_error'),
                    'ci_lower': res.get('ci_lower'),
                    'ci_upper': res.get('ci_upper'),
                    'total_time': res.get('total_time'),
                    'status': 'success'
                }
                writer.writerow(row)
                csvfile.flush()
                
            except Exception as e:
                print(f"  FAILED {method_name}: {e}")
                traceback.print_exc()
                row = {
                    'model': model_name,
                    'method': method_name,
                    'status': f'error: {str(e)}',
                    'probability': -1, 'std_error': -1, 'ci_lower': -1, 'ci_upper': -1, 'total_time': -1
                }
                writer.writerow(row)
                    
    print("\nBenchmark Finished. Results saved to", output_file)

if __name__ == "__main__":
    main()

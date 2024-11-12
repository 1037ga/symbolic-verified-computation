import subprocess

# 実行するスクリプトのファイル名
script_name = 'interval_n_vars_test.py'

# 渡したい入力のリスト
inputs1 = [i for i in range(1,7)]
inputs2 = [i for i in range(1,20)]

for input_value1 in inputs1:
    for input_value2 in inputs2:
        # コマンドを構築
        command = f'python3 {script_name} {input_value1} {input_value2}'
        print(command)
        # サブプロセスを起動してスクリプトを実行
        subprocess.run(command, shell=True)

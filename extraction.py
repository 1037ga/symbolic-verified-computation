def extraction(filename):
    detas = []
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if i % 9 == 7:
                # deta = (line[len('       '):].split('+')[1])
                deta = line
                deta = deta[5:len(deta)-1]
                detas.append(deta)
    return detas

if __name__ == "__main__":
    # ファイル名
    filename = 'output7.txt'
    # 名前を抽出
    extracted_detas = extraction(filename)
    # 結果を出力
    for deta in extracted_detas:
        print(deta)

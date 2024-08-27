import os


def main():
    o = open("/home/avraham/MaruvkaLab/Texas/efficient_run/o.full.mut.tsv", 'r')
    n = open("/home/avraham/MaruvkaLab/Texas/efficient_run/n.full.mut.tsv", 'r')
    o.readline()
    n.readline()
    i=0
    while True:
        i+=1
        o_line_tokens = o.readline().split("\t")
        if o_line_tokens[0] == '':
            break
        n_line_tokens = n.readline().split("\t")
        for o_t, n_t in zip(o_line_tokens, n_line_tokens):
            if o_t!=n_t:
                try:
                    o_tf = float(o_t)
                    n_tf = float(n_t)
                    if abs(o_tf-n_tf) > 0.01:
                        print(f"Line: {i}\nOld: {o_t}\nNew: {n_t}\n***************************************")
                except:
                    print(f"Line: {i}\nOld: {o_t}\nNew: {n_t}\n***************************************")

    o.close()
    n.close()



if __name__ == '__main__':
    main()
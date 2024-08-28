

def create_tsv(header, oline, nline):
    with open("/home/avraham/MaruvkaLab/Texas/efficient_run/full/diffline.tsv", 'w+') as diffline:
        diffline.write(header)
        diffline.write(oline)
        diffline.write(nline)


def main():
    o = open("/home/avraham/MaruvkaLab/Texas/efficient_run/full/ofull.full.mut.tsv", 'r')
    n = open("/home/avraham/MaruvkaLab/Texas/efficient_run/full/nfull.full.mut.tsv", 'r')
    oline = o.readline()
    header=oline
    nline = n.readline()
    i=0
    while oline!="":
        i+=1
        otokens = oline.split('\t')
        ntokens = nline.split('\t')
        ocall = otokens[48]
        ncall = ntokens[48]
        if otokens[48]!=ntokens[48]:
            print("Line_num" + str(i))
            print(otokens[1])
            print(ocall)
            print(ncall)
            croc=1
        oline = o.readline()
        nline = n.readline()


main()
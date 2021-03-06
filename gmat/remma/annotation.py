import re
import sys


def gtf_to_gene_info(gtf_file):
    fout = open(gtf_file + '.gene_info', 'w')
    fin = open(gtf_file)
    for line in fin:
        if '#' in line:
            continue
        arr = line.split()
        if arr[2] == 'gene':
            searchObj = re.search('gene_id\s+"(.+?)".+gene_name\s+"(.+?)"', line, re.I)
            gene_id = searchObj.group(1)
            gene_name = searchObj.group(2)
            gene_info = ' '.join([arr[0], arr[3], arr[4], arr[6], gene_id, gene_name])
            fout.write(gene_info + '\n')
    fin.close()
    fout.close()


def annotation_snp_pos(res_file, bed_file, p_cut=1, dis=0, ld_file=None, r2=0.2):
    """
    add the snp information for the epistasis result file; select top significant SNP pairs based on p cut value;
    select SNP pairs bases on their distance.
    :param res_file: The result file for epistasis file. The first two columns are the orders for the SNP pairs, and
    the last column is the p values.
    :param bed_file: the prefix for the plink binary file
    :param p_cut: p cut value, default value is 1.
    :param dis: the min distance between SNP pairs, the default value is -1.
    :return: 0
    """
    snp_info = {}
    order = -1
    with open(bed_file + '.bim', 'r') as fin:
        for line in fin:
            order += 1
            arr = line.split()
            snp_info[str(order)] = ' '.join(arr)
    with open(res_file, 'r') as fin, open(res_file + '.anno', 'w') as fout:
        line = fin.readline()
        arr = line.split()
        fout.write(' '.join([arr[0], 'snp0_chr', 'snp0_ID', 'snp0_cm', 'snp0_bp', 'snp0_allele1', 'snp0_allele2',
                         arr[1], 'snp1_chr', 'snp1_ID', 'snp1_cm', 'snp1_bp', 'snp1_allele1', 'snp1_allele2']))
        fout.write(' ')
        fout.write(' '.join(arr[2:]))
        fout.write('\n')
        for line in fin:
            arr = line.split()
            snp0 = snp_info[arr[0]].split()
            snp1 = snp_info[arr[1]].split()
            if float(arr[-1]) <= p_cut and (snp0[0] != snp1[0] or abs(float(snp0[3]) - float(snp1[3])) > dis):
                fout.write(' '.join([arr[0], snp_info[arr[0]], arr[1], snp_info[arr[1]]]))
                fout.write(' ')
                fout.write(' '.join(arr[2:]))
                fout.write('\n')
    if ld_file is not None:
        ld_id= {}
        with open(ld_file, 'r') as fin:
            fin.readline()
            for line in fin:
                arr = line.split()
                if float(arr[-1]) > r2:
                    ld_id[' '.join([arr[2], arr[5]])] = 1
                    ld_id[' '.join([arr[5], arr[2]])] = 1
        with open(res_file + '.anno', 'r') as fin, open(res_file + '.anno.ld', 'w') as fout:
            line = fin.readline()
            fout.write(line)
            for line in fin:
                arr = line.split()
                snp = ' '.join([arr[2], arr[9]])
                if snp not in ld_id:
                    fout.write(line)
    return 0


def annotation_snp_nearest_gene(bed_file, gene_file, max_distance=150000):
    """
    ??????SNP????????????
    :param bed_file: plink??????
    :param gene_file: ???????????????????????????????????????????????????????????????????????????
    :param max_distance: SNP???????????????????????????
    :return:
    """
    gene_info = {}
    fin = open(gene_file)
    for line in fin:
        arr = line.split()
        try:
            gene_info[arr[0]].append(arr[:])
        except Exception as e:
            del e
            gene_info[arr[0]] = [arr[:]]
    fin.close()
    fout = open(bed_file + '.nearby_genes', 'w')
    fin = open(bed_file + '.bim')
    for line in fin:
        snp_info = line.strip()
        arr = line.split()
        snp_pos = int(arr[3])
        for genei in gene_info[arr[0]]:
            start = int(genei[1])
            end = int(genei[2])
            dist1 = snp_pos - start
            dist2 = snp_pos - end
            if dist1 > 0 and dist2 < 0:
                fout.write(snp_info + ' ' + ' '.join(genei) + ' within' + '\n')
            else:
                distance = min(abs(dist1), abs(dist2))
                if distance < max_distance:
                    fout.write(snp_info + ' ' + ' '.join(genei) + ' ' + str(distance) + '\n')
    fin.close()
    fout.close()



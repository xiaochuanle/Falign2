import os
import shutil
import sys

add_error_script_path = sys.argv[1]
ds_script_path = sys.argv[2]
input_reads_path = sys.argv[3]
output_dir = sys.argv[4]
output_fastq_path = sys.argv[5]

FQDIR=f"{output_dir}/sim-fqs"
if os.path.exists(FQDIR):
    shutil.rmtree(FQDIR)
os.mkdir(FQDIR)

infp = open(input_reads_path, 'r')

readBatchSize = 1000

def make_fa_read_path(id):
    return f"{output_dir}/{id}.fasta"

def make_fq_read_path(id):
    return f"{output_dir}/sim-fqs/{id}.fastq"

def make_next_read_batch():
    id = 0
    while True:
        hdr = infp.readline()
        if not hdr:
            break
        hdr = hdr.strip()
        assert hdr[0] == '>'
        
        seq = infp.readline()
        assert seq is not None
        seq = seq.strip()

        output_path = make_fa_read_path(id)
        out = open(output_path, 'w')
        print(hdr, file = out)
        print(seq, file = out)
        out.close()

        id += 1
        if id == readBatchSize:
            break 
    return id 

def main():
    out = open(output_fastq_path, 'w')
    while True:
        numReads = make_next_read_batch()
        if numReads == 0:
            break
        print(f"Dump {numReads} reads")

        cmd = f"{add_error_script_path} {ds_script_path} {output_dir} {FQDIR} {numReads}"
        r =os.system(cmd)
        r = r>>8
        if r != 0:
            sys.exit(1)

        for i in range(numReads):
            fafn = make_fa_read_path(i)
            fafp = open(fafn, 'r')
            fahdr = fafp.readline()
            fafp.close()
            assert fahdr is not None
            fahdr = fahdr.strip()
            fahdr = fahdr[1:]

            fqfn = make_fq_read_path(i)
            fqfp = open(fqfn, 'r')
            hdr = fqfp.readline()
            hdr = hdr.strip()
            hdr = hdr.split(' ')[0]
            seq = fqfp.readline()
            seq = seq.strip()
            plus = fqfp.readline()
            plus = plus.strip()
            qual = fqfp.readline()
            qual = qual.strip()
            fqfp.close()

            newhdr = hdr + ' ' + fahdr
            print(newhdr, file = out)
            print(seq, file = out)
            print(plus, file = out)
            print(qual, file =out)

        
        #break

    out.close()

if __name__ == '__main__':
    main()

#!/usr/bin/python3

import sys

from flytekit import kwtypes, task, workflow
from flytekit.extras.tasks.shell import ShellTask

t1 = ShellTask(
    name="task_1",
    debug=True,
    script="""
    #!/usr/bin/bash

    mkdir results
    mkdir results/fastqc_reports
    fastqc {inputs.srr1} {inputs.srr2} --outdir=results/fastqc_reports/
    mkdir temp
    echo "BWA INDEX\n"
    bwa index {inputs.ref}
    echo "FINISHED BWA INDEX\n"
    echo "BWA MEM\n"
    bwa mem {inputs.ref} {inputs.srr1} {inputs.srr2} -o sample.sam
    echo "FINISHED BWA MEM\n"
    echo "SAMTOOLS VIEW\n"
    samtools view -b -o sample.bam sample.sam
    echo "FINISHED SAMTOOLS VIEW\n"
    echo "SAMTOOLS FLAGSTAT\n"
    samtools flagstat sample.bam > results/flagstat_result.txt
    """,
    inputs=kwtypes(srr1=str, srr2=str, ref=str)
)

t_ok = ShellTask(
    name="task_1",
    debug=True,
    script="""
    #!/usr/bin/bash
    
    echo "OK"
    echo "SAMTOOLS SORT AND INDEX"
    samtools sort -o sorted.bam sample.bam
    samtools index sorted.bam
    echo "FINISHED SAMTOOLS SORT AND INDEX"
    echo "FREEBAYES"
    freebayes -f reference.fna -b sorted.bam --vcf results/sample.vcf
    echo "NOT OK"
    """
)

t_not_ok = ShellTask(
    name="task_1",
    debug=True,
    script="""
    #!/usr/bin/bash

    echo "NOT OK"
    """
)


def get_num():
    with open("results/flagstat_result.txt") as f:
        ss = f.read().split('\n')
    n1 = ss[7].find('(') + 1
    n2 = ss[7].find('%')
    res = float(ss[7][n1:n2])
    return res


@task
def run_first_part() -> None:
    t1(srr1=sys.argv[1], srr2=sys.argv[2], ref=sys.argv[3])


@task
def run_second_part() -> None:
    evaluation = get_num() > 90
    if evaluation:
        t_ok()
    else:
        t_not_ok()


@workflow
def wf() -> None:
    run_first_part()
    run_second_part()


if __name__ == "__main__":
    wf()

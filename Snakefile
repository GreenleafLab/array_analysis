rule CPfluor2CPseries:
    input:
        {condition}.CPfluor
        NNNlib2.map
    output:
        {condition}.CPseries.csv

rule 
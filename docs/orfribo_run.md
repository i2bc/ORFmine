## Launch the analysis

Be sure that your configuration.yaml file is [completed](./orfribo_configuration.md) and stored in the /workdir/ directory of the container. Then, from the /workdir/ directory, please type the following command:

``` bash
orfribo CPU_NUMBER MEMORY_AMOUNT
```
Where *CPU_NUM* is the number of CPU/threads and MEMORY_AMOUNT is the
amount of memory (in Giga-Bytes) that can be used for the analysis.

Example for 6 CPUs and 12GB memory use:
``` bash
orfribo 6 12
```

One step of the analysis is done by the [riboWaltz tool](https://github.com/LabTranslationalArchitectomics/riboWaltz) which can need high RAM resources depending on your data. If the pipeline crashes without any obvious reason, it might be because of a lack of memory available.

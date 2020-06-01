# Retrieving TARGET-OS data from GDC Commons

This is the procedure I used to get the data from GDC Commons:

Visit the (GDC site for the TARGET-OS RNAseq data)[https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TARGET-OS%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_category%22%2C%22value%22%3A%5B%22Transcriptome%20Profiling%22%5D%7D%7D%5D%7D&searchTableTab=cases].

On the 'Files' tab on the left, select the HTSeq-counts check box under the 'Workflow Type'. Add add of the associated files to your Cart using the 'Add all files to cart' button. Once the files are selected, click the 'Download manifest' button and save the file to the location where you will store the files. After downloading the manifest file, click the 'Biospecimen' and 'Clinical' files to download the associated metadata files.

Now we need to get the GDC transfer tool. The tool can be downloaded from (here)[https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/]. You can save it wherever you would like, but I just dropped in the same location as the manifest file for simplicity.

I am on a Windows 10 system, so from here I opened command prompt, navigated to the directory where I stored the files, and downloaded the files in the manifest file.

```
>cd OneDrive\Documents\bccrc\projectsRepository\sorensenLab\relatedToYbx1\20200601_explorationYbx1SurvivalOsTarget
>gdc-client.exe download -m gdc_manifest.2020-05-29.txt
```

The files are each in their own folder as a .gz archive. To extract them all at once, I used my Ubuntu subsytem terminal on my Windows 10 machine.

```
$cd /mnt/c/Users/chris/OneDrive/Documents/bccrc/projectsRepository/sorensenLab/relatedToYbx1/20200601_explorationYbx1SurvivalOsTarget/
$find . -name '*.gz' | while read filename; do gzip -dv $filename; done;
```

You should now have a folder containing text files for counts, per patient, in the original directory where everything was saved. We will now process the count data in R.

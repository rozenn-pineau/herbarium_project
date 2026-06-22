### Doanloading the data from Novogene usinf lftp
```
lftp -c 'set sftp:auto-confirm yes;set net:max-retries 20;open sftp://X202SC26058003-Z01-F002:uD2dCtY9@usftp22.novogene.com; mirror --verbose --use-pget-n=8 -c'
```
### Checking the download
```
md5sum -c MD5.txt
```

Every file looked good.

### Demultiplexing the fastq

(a) Match barcode to sequence:

I used a custom R script to match the barcodes i5 and i7 to their corresponding sequence (match_barcode_sequence.R).

(b) demultiplex:

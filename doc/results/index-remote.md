 Results for RSQF
 ```
Running file_fpr_test for data/chicagoB64.txt (1700741 unique lines) with A/S=1000, load=0.95, rem_size=8...
Initializing filter and pulling lines from file...
Inserting 1945 items into filter...
Querying the filter with the remaining lines...
False positives: 63482 (0.1399575312383944%), repeated: 61210 (0.13494849700863432%)
False negatives: 0 (0%)
Running file_fpr_test for data/chicagoB64.txt (1700741 unique lines) with A/S=10000, load=0.95, rem_size=8...
Initializing filter and pulling lines from file...
Inserting 243 items into filter...
Querying the filter with the remaining lines...
False positives: 84295 (0.18583657444121104%), repeated: 81145 (0.17889209126320743%)
False negatives: 0 (0%)

 ```

Results for index-remote extension-based AQF
```
Running file_fpr_test for data/chicagoB64.txt (1700741 unique lines) with A/S=1000, load=0.95, rem_size=8...
Initializing filter and pulling lines from file...
Inserting 1945 items into filter...
Querying the filter with the remaining lines...
False positives: 3969 (0.008750377138168102%), repeated: 1860 (0.004100705839504326%)
False negatives: 0 (0%)
Running file_fpr_test for data/chicagoB64.txt (1700741 unique lines) with A/S=10000, load=0.95, rem_size=8...
Initializing filter and pulling lines from file...
Inserting 243 items into filter...
Querying the filter with the remaining lines...
False positives: 17880 (0.039418209277049095%), repeated: 14838 (0.032711822665148466%)
False negatives: 0 (0%)
```

Results for selector-based AQF (C code)
```
Running AQF integration tests...
Running asymmetric test with ../data/chicagoB64.txt (1700741 uniques), using first 1945 unique elements as set...
Inserting 1945 members into filter...done.
Running queries...done.
False negatives: 0
False positives: 14407; Repeated: 7175

Total Time: 12.254249
Insert: 0.002564
Query: 12.250586
0.000001	 0.000000

Running asymmetric test with ../data/chicagoB64.txt (1700741 uniques), using first 243 unique elements as set...
Inserting 243 members into filter...done.
Running queries...done.
False negatives: 0
False positives: 33595; Repeated: 26076

Total Time: 10.157856
Insert: 0.000308
Query: 10.157445
0.000001	 0.000000

Running asymmetric test with ../data/sanjose64.txt (2193052 uniques), using first 31129 unique elements as set...
Inserting 31129 members into filter...done.
Running queries...done.
False negatives: 0
False positives: 8626; Repeated: 383

Total Time: 14.345554
Insert: 0.056019
Query: 14.270929
0.000002	 0.000000
```

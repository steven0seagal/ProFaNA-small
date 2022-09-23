# ProFaNA-small
## Changes
- connection to postgresql
- smaller matrix creating and 6list reading

## Benchmarks for gene list as an input (7426 genes)
1. noSQL (total_time -> 00:08:09)
 
   - DLOK TIME START 2022-09-13 09:37:11.581376
   - DLOK TIME OVER 2022-09-13 09:44:05.887296
   - MATRIX CREATING 2022-09-13 09:44:06.304499
   - MATRIX OVER 2022-09-13 09:45:07.080097
   - Wilcoxon start 2022-09-13 09:45:07.080122
   - Wilcoxon OVER 2022-09-13 09:45:15.348673
   - Collecting data START 2022-09-13 09:45:15.348723
   - Collecting data OVER 2022-09-13 09:45:17.596108
   - Customizing_data START 2022-09-13 09:45:17.596143
   - Customizing_data OVER 2022-09-13 09:45:20.179166


2. SQL at genome findings (total_time -> 00:03:07 )

   - DLOK TIME START 2022-09-13 12:34:00.454407
   - DLOK TIME OVER 2022-09-13 12:35:57.087847
   - MATRIX CREATING 2022-09-13 12:35:57.505363
   - MATRIX OVER 2022-09-13 12:36:55.579488
   - Wilcoxon start 2022-09-13 12:36:55.579515
   - Wilcoxon OVER 2022-09-13 12:37:03.794699
   - Collecting data START 2022-09-13 12:37:03.794724
   - Collecting data OVER 2022-09-13 12:37:04.807227
   - Customizing_data START 2022-09-13 12:37:04.807252
   - Customizing_data OVER 2022-09-13 12:37:07.115529


3. SQL at genome finding and genome size searching ( total time -> 00:03:40)

   - OK TIME START 2022-09-13 13:29:25.616570
   - DLOK TIME OVER 2022-09-13 13:31:57.303319
   - MATRIX OVER 2022-09-13 13:32:58.041552 
   - Wilcoxon start 2022-09-13 13:32:58.041578 
   - Wilcoxon OVER 2022-09-13 13:33:06.109163 
   - Collecting data START 2022-09-13 13:33:06.109186 
   - Collecting data OVER 2022-09-13 13:33:07.105116
   - Customizing_data START 2022-09-13 13:33:07.105142 
   - Customizing_data OVER 2022-09-13 13:33:09.382293


4. SQLs and smaller matrix creating (total time -> 00:02:57)

   - DLOK TIME START 2022-09-13 14:39:25.768522
   - DLOK TIME OVER 2022-09-13 14:42:09.930909
   - MATRIX CREATING 2022-09-13 14:42:10.368518
   - MATRIX OVER 2022-09-13 14:42:10.959046
   - Wilcoxon start 2022-09-13 14:42:10.959072
   - Wilcoxon OVER 2022-09-13 14:42:18.872992
   - Collecting data START 2022-09-13 14:42:18.873017
   - Collecting data OVER 2022-09-13 14:42:19.847282
   - Customizing_data START 2022-09-13 14:42:19.847310
   - Customizing_data OVER 2022-09-13 14:42:22.099046

5. THE BEST -- One SQL in genome findings and smaller matrix creating (total time --> 00:02:28)

   - DLOK TIME START 2022-09-13 14:53:10.276174
   - DLOK TIME OVER 2022-09-13 14:55:25.349914
   - MATRIX CREATING 2022-09-13 14:55:25.775047
   - MATRIX OVER 2022-09-13 14:55:26.371425 
   - Wilcoxon start 2022-09-13 14:55:26.371451 
   - Wilcoxon OVER 2022-09-13 14:55:34.846980
   - Collecting data START 2022-09-13 14:55:34.847005
   - Collecting data OVER 2022-09-13 14:55:35.923965
   - Customizing_data START 2022-09-13 14:55:35.923990
   - Customizing_data OVER 2022-09-13 14:55:38.563050

6. Same as above but with indexed database ( total time 00:02:08)

   - DLOK TIME START 2022-09-22 21:17:38.767107
   - DLOK TIME OVER 2022-09-22 21:19:33.845042
   - MATRIX CREATING 2022-09-22 21:19:34.262756
   - MATRIX OVER 2022-09-22 21:19:34.840180
   - Wilcoxon start 2022-09-22 21:19:34.840206
   - Wilcoxon OVER 2022-09-22 21:19:43.064248
   - Collecting data START 2022-09-22 21:19:43.064270
   - Collecting data OVER 2022-09-22 21:19:44.063362
   - Customizing_data START 2022-09-22 21:19:44.063386
   - Customizing_data OVER 2022-09-22 21:19:46.572826

## Benchmark for domain(pfam02696) as an input with 1000 bp neighborhod on all_database

1. No changes (total time --> 00:27:12)

    - DLOK TIME START 2022-09-13 15:58:18.492637
    - DLOK TIME OVER 2022-09-13 16:15:02.575725 
    - MATRIX CREATING 2022-09-13 16:15:03.102840 
    - MATRIX OVER 2022-09-13 16:25:05.018804
    - Wilcoxon start 2022-09-13 16:25:05.018834
    - Wilcoxon OVER 2022-09-13 16:25:21.528457
    - Collecting data START 2022-09-13 16:25:21.528483
    - Collecting data OVER 2022-09-13 16:25:26.701745
    - Customizing_data START 2022-09-13 16:25:26.701778
    - Customizing_data OVER 2022-09-13 16:25:30.297668
2. Grep + use gene_list approach ( total time -->00:13:18)

    - START 2022-09-13 20:37:15.081420
    - Checkpoint #1 Data loaded 2022-09-13 20:37:54.781407
    - DLOK TIME START 2022-09-13 20:37:54.790883
    - DLOK TIME OVER 2022-09-13 20:50:03.023851
    - MATRIX CREATING 2022-09-13 20:50:03.581646
    - MATRIX OVER 2022-09-13 20:50:07.707844
    - Wilcoxon start 2022-09-13 20:50:07.707870
    - Wilcoxon OVER 2022-09-13 20:50:22.816301
    - Collecting data START 2022-09-13 20:50:22.816330
    - Collecting data OVER 2022-09-13 20:50:30.763353
    - Customizing_data START 2022-09-13 20:50:30.763385
    - Customizing_data OVER 2022-09-13 20:50:33.816410
    - END 2022-09-13 20:50:33.857694

3. Same as above but with indexex (total time -->00:11:19)

   - START 2022-09-22 21:23:54.276017
   - Checkpoint #1 Data loaded 2022-09-22 21:24:34.011121
   - DLOK TIME START 2022-09-22 21:24:34.020479
   - DLOK TIME OVER 2022-09-22 21:34:44.364148
   - MATRIX CREATING 2022-09-22 21:34:44.927492
   - MATRIX OVER 2022-09-22 21:34:49.145238
   - Wilcoxon start 2022-09-22 21:34:49.145259
   - Wilcoxon OVER 2022-09-22 21:35:03.876001
   - Collecting data START 2022-09-22 21:35:03.876026
   - Collecting data OVER 2022-09-22 21:35:12.039912
   - Customizing_data START 2022-09-22 21:35:12.039946
   - Customizing_data OVER 2022-09-22 21:35:15.170198
   - END 2022-09-22 21:35:15.209903

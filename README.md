# DataScience Master Thesis

This project aims to compare, in therms of quality and performance, three of the most used Samtools commands with the result of correspondant SQL queries.
Samtools commands chosen for the comparison are the following: consensus, mpileup, coverage

The file adopted in the project ENCFF469GFN.bam is freely dowloadable from its source: https://www.encodeproject.org/files/ENCFF469GFN/

The project is structured by the following steps:

1. PostrgeSQL schema implementation                           --> DBcreator.py
2. BAM file insertion into the database schema                --> DBinsertion.py
3. SQL queries creation                                       --> consensus.sql, coverage.sql, mpileup.sql
4. SQL index creation                                         --> index_creation.sql
5. CSV files have been produced for performance evaluation    --> DBinsertion_performance.csv, query_performance_log.csv, samtools_performance_log.csv
6. CSV files have been processed for statistical evaluations  --> processing_DBinsertion_performance.ipynb, processing_query_performance_log.ipynb, processing_samtools_performance_log.ipynb

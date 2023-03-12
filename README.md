This project aims to compare, in therms of quality and performance, three of the most used SAMtools commands with the result given by respective sperimental programs.
Samtools commands chosen for the comparison are the following: consensus, mpileup, coverage.

3 versions of the programs have been developed: 2 versions based on SQLite sperimenting 2 different schema versions, 1 version built directly on PostgreSQL.
The performance comparison will be studied between the sperimental version and Samtools.
Imput datasets are SAM files.

The project is structured by the following steps:

The directory is structured as follows:
- Perfomance: here are found the scripts used to evaluate the perfomrance of the different versions in therms of quality, speed and RAM consumption.
- V_PostgreSQL: commands scripts in pl/pgSQL, database creation and population from SAM file.
- V1_SQlite: python programs corresponding to SAMtools commands and database creation and population from SAM file. Basic database schema is adopted.
- V2_SQLite: python programs corresponding to SAMtools commands and database creation and population from SAM file. Alternative database schema is adopted.

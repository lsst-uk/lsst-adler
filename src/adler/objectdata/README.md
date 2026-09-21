adler_data_schema.py

This schema-parameter table is used to map the parameter name in a given schema to the name that adler expects.

Column: {schema} - names of the parameter in {schema}
Column: {schema}_table - the table name that each parameter is found

Available {schema}:
adler - parameters refer to adler class attributes, table to their adler class.
dp03_catalogs_10yr - access via RSP
dp1 - access via RSP
dp2 - access via RSP
MPC - access via https://b612.ai/rubin-mpc-downloads/ (missing additional ephemerides) or via SSSC release on RSP /home/mschwamb/RubinSSPviaMPC/

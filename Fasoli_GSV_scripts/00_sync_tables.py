#!/usr/bin/env python

# from multiprocessing import Pool

from google.cloud import bigquery


client = bigquery.Client()


dataset_ref = "gsv-aq-x.saltlakecity"
tables = [
    # "2_2b_205",
    # "2_gps",
    # "2_lgr_mgga",
    # "2_licor_7000",
    # "2_magee_ae33",
    # "2_magic_200p",
    # "2_picarro_g2401",
    # "2_teledyne_t200",
    # "2_teledyne_t500u",
    # "2_thermo_pdr_1500",
    # "2_vaisala_wxt536",
    # "3_primary_keys",
    # "3_primary_keys_nulled",
    "4_primary_keys_uuids"
]


def fetch_table(table: str):
    query = f"SELECT * FROM `{dataset_ref}.{table}`"
    df = client.query(query).result().to_dataframe()
    df = df.sort_values(by=["system_id", "time"]).reset_index()
    # df.to_csv(f"data/tables/{table}.csv", index=False)
    df.to_feather(f"data/tables/{table}.feather")
    # df.to_pickle(f"data/tables/{table}.p")
    # df.to_parquet(f"data/tables/{table}.parquet")
    return table


if __name__ == "__main__":
    # pool = Pool(len(tables))
    # pool.map(fetch_table, tables)

    for table in tables:
        fetch_table(table)

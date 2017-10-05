CREATE TABLE species
(
    taxid INTEGER primary key,
    domain TEXT,
    organism_name TEXT,
    ftp_path TEXT,
    ftp_basename TEXT,
    genetic_code INTEGER
);

CREATE TABLE flow
(
    taxid INTEGER primary key,
    download INTEGER DEFAULT 0,
    preprocess INTEGER DEFAULT 0,
    shuffle INTEGER DEFAULT 0,
    coverage INTEGER DEFAULT 0
);

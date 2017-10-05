import sqlite3


class DbController:
    def __init__(self, dbFilepath):

        def dict_factory(cursor, row):
            d = {}
            for idx, col in enumerate(cursor.description):
                d[col[0]] = row[idx]
            return d

        self.dbFilepath = dbFilepath
        self.con = sqlite3.connect(self.dbFilepath)
        self.con.row_factory = dict_factory
        self.cur = self.con.cursor()

    def __del__(self):
        self.con.close()

    def execute(self, query, arg=None):
        try:
            if arg is None:
                self.cur.execute(query)
            else:
                self.cur.execute(query, arg)
            self.con.commit()
        except (sqlite3.IntegrityError, sqlite3.OperationalError) as e:
            print(e)
            return False
        return True

    def insert_species(self, taxid, domain, organism_name, ftp_path, ftp_basename, genetic_code):
        query = 'INSERT INTO species VALUES(?, ?, ?, ?, ?, ?)'
        arg = (taxid,  domain, organism_name, ftp_path, ftp_basename, genetic_code)
        success = self.execute(query, arg)

        if success:
            query = 'INSERT OR REPLACE INTO flow(taxid) VALUES(?)'
            arg = (taxid,)
            success = self.execute(query, arg)
            return success
        else:
            return False

    def get_target(self, process):
        assert process in ["download", "preprocess", "shuffle", "coverage"]

        query = 'SELECT taxid FROM flow WHERE {}=0'.format(process)
        success = self.execute(query)
        if success:
            ret = [_["taxid"] for _ in self.cur.fetchall()]
            return ret
        else:
            return []

    def get_row(self, taxid):
        query = 'SELECT * FROM species WHERE taxid=?'
        arg = (taxid, )
        success = self.execute(query, arg)

        if success:
            return self.cur.fetchone()
        else:
            return None

    def mark_as_done(self, taxid, process):
        assert process in ["download", "preprocess", "shuffle", "coverage"]
        query = 'UPDATE flow SET {}=1 WHERE taxid=?'.format(process)
        arg = (taxid, )
        success = self.execute(query, arg)
        return success

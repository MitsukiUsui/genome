from myutil.myutil import DbController

class AltorfDbController(DbController):
    def __init__(self, dbFilepath):
        super().__init__(dbFilepath)
    
    def insert_species(self, taxid, domain, organism_name, ftp_path, ftp_basename, genetic_code):
        """
        Try to insert new row into species table. If succeeded, initialize flow table.
        """
        
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

    def get_target(self, process=None):
        """
        given process name, return a list of target taxid
        """
        
        if process is None:
            query = 'SELECT taxid FROM flow'
        else:
            assert process in self.__columns("flow")
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

    def mark_as_done(self, taxid_lst, process):
        assert process in ["download", "preprocess", "shuffle", "coverage"]
        query = 'UPDATE flow SET {}=1 WHERE taxid in ({})'.format(process, ','.join([str(taxid) for taxid in taxid_lst]))
        success = self.execute(query)
        return success
    
if __name__=='__main__':
    adc = AltorfDbController("../db/altorf.db")
    print(adc.get_target()[:5])

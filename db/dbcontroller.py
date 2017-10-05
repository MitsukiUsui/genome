import sqlite3

class DbController:
    def __init__(self, dbFilepath):
        self.dbFilepath = dbFilepath
        self.con=sqlite3.connect(self.dbFilepath)
        self.cur=self.con.cursor()
    
    def __del__(self):
        self.con.close()
            
    def execute(self, query, arg):
        try:
            self.cur.execute(query, arg)
            self.con.commit()
        except (sqlite3.IntegrityError, sqlite3.OperationalError) as e:
            print (e)
            return False
        return True
        
    def insert_species(self, taxid, domain, organism_name, ftp_path, ftp_basename, genetic_code):        
        query='INSERT INTO species VALUES(?, ?, ?, ?, ?, ?)'
        arg=(taxid,  domain, organism_name, ftp_path, ftp_basename, genetic_code)
        success=self.execute(query, arg)
        
        if success:
            query='INSERT OR REPLACE INTO flow(taxid) VALUES(?)'
            arg=(taxid,)
            success = self.execute(query, arg)
            return success
        else:
            return False

import sqlite3

class StringIndex:
    _createTable = """
        create table stringindex (
            id integer primary key asc autoincrement,
            str text,
            cnt integer default 1
         );
    """
    _index1 = """
        create unique index index1 on stringindex (str)
        """
    _insert = """
        insert into stringindex (id,str,cnt) values (NULL,?,?)
    """
    _find = """
        select id from stringindex where str = ?
    """
    _match = """
        select id from stringindex where str regexp ?
    """
    _increment = """
        update stringindex set cnt = cnt + 1 where id = ?
    """
    _strings = """
        select str from stringindex where cnt > 0 order by str
    """
    _ids = """
        select id from stringindex where cnt > 0 order by id
    """
    _byid = """
        select id,str,cnt from stringindex where cnt > 0 order by id
    """
    _bystr = """
        select id,str,cnt from stringindex where cnt > 0 order by str
    """
    _bycnt = """
        select id,str,cnt from stringindex where cnt > 0 order by cnt
    """
    _counts = """
        select id,cnt from stringindex where cnt > 0 order by id
    """
    _count = """
        select cnt from stringindex where id = ?
    """
    _string = """
        select str from stringindex where id = ?
    """
    _len = """
        select count(*) from stringindex where cnt > 0;
    """
    def __init__(self):
        conn = sqlite3.connect(':memory:',isolation_level='DEFERRED')
        conn.text_factory = str
        conn.execute(self._createTable)
        conn.execute(self._index1)
        self.conn = conn
    def id(self,s):
        for r in self.conn.execute(self._find,(str(s),)):
            return r[0]
        return None
    def has(self,s):
        return (self.id(s) != None)
    def get(self,s):
        try:
            cursor = self.conn.cursor()
            cursor.execute(self._insert,(str(s),0))
            id = cursor.lastrowid
        except sqlite3.IntegrityError:
            id = self.id(s)
        return id
    def add(self,s):
        try:
            cursor = self.conn.cursor()
            cursor.execute(self._insert,(str(s),1))
            id = cursor.lastrowid
        except sqlite3.IntegrityError:
            id = self.id(s)
            self.conn.execute(self._increment,(id,))
        return id
    def match(self,s):
        for r in self.conn.execute(self._match,(str(s),)):
            yield r[0]
    def update(self,l):
        for s in l:
            self.add(s)
    def strings(self):
        for r in self.conn.execute(self._strings):
            yield r[0]
    def ids(self):
        for r in self.conn.execute(self._ids):
            yield r[0]
    def count(self,id):
        for r in self.conn.execute(self._count,(id,)):
            return r[0]
    def string(self,id):
        for r in self.conn.execute(self._string,(id,)):
            return r[0]
    def counts(self):
        for r in self.conn.execute(self._counts):
            yield r[0],r[1]
    def byid(self):
        for r in self.conn.execute(self._byid):
            yield r
    def bystring(self):
        for r in self.conn.execute(self._bystr):
            yield r
    def bycount(self):
        for r in self.conn.execute(self._bycount):
            yield r
    def __iter__(self):
        return self.ids()
    def index(self,s):
        return self.add(s)
    def __len__(self):
        for r in self.conn.execute(self._len):
            return r[0]
    def __str__(self):
        lines = []
        for i,s,c in self.byid():
            lines.append('\t'.join(map(str,[i,s,c])))
        return '\n'.join(lines)

class BipartiteStringIndex:
    _createTable = """
        create table edges (a int, b int, cnt int default 1);
    """
    _index1 = """
        create unique index index1 on edges (a,b);
    """
    _index2 = """
        create unique index index2 on edges (b,a);
    """
    _insert = """
        insert into edges (a,b) values (?,?)
    """
    _increment = """
        update edges set cnt = cnt + 1 where a=? and b=?
    """
    _out = """
        select b from edges where a = ?
    """
    _into = """
        select a from edges where b = ?
    """
    def __init__(self,indexA,indexB):
        self.a = indexA
        self.b = indexB
        conn = sqlite3.connect(':memory:',isolation_level='DEFERRED')
        conn.execute(self._createTable)
        conn.execute(self._index1)
        conn.execute(self._index2)
        self.conn = conn

    def add(self,u,v):
        aid = self.a.index(u)
        bid = self.b.index(v)
        try:
            self.conn.execute(self._insert,(aid,bid))
        except ValueError:
            self.conn.execute(self._increment,(aid,bid))

    def out(self,u):
        return set(r[0] for r in self.conn.execute(self._out,u))

    def into(self,v):
        return set(r[0] for r in self.conn.execute(self._into,v))

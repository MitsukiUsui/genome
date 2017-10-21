class Phase:
    def __add(self,p):
        i = int(p) + 1
        if p[0]=='+':
            return "+{}".format(i % 3)
        elif p[0]=='-':
            if i == 1:
                i -= 3
            return "-{}".format(abs(i))
    
    def __rev(self, p):
        if p[0]=='+':
            return "-{}".format(abs(int(p)))
        elif p[0]=='-':
            return "+{}".format(abs(int(p)))

    
    def __operate(self, p1, p2):
        """
        p1に対してp2の関係にあるphaseを返す
        """
        for op in self.ph2op[p2]:
            if op=='a':
                p1=self.__add(p1)
            elif op=='r':
                p1=self.__rev(p1)
        return p1
    
    def __init__(self):       
        #define phase and its operation notation.
        self.phase_lst=["+0", "+1", "+2", "-0", "-1", "-2"]
        self.ph2op={"+0":[], "+1":['a'], "+2":['a', 'a'],
                    "-0":['r'], "-1":['a', 'r'], "-2":['a', 'a', 'r']}
        
        self.operate_dct={}
        for p1 in self.phase_lst:
            for p2 in self.phase_lst:
                self.operate_dct[(p1, p2)]=self.__operate(p1, p2)
                
        self.relative_dct={}
        for p1 in self.phase_lst:
            for p2 in self.phase_lst:
                p3=self.operate_dct[(p1, p2)]
                self.relative_dct[(p1, p3)]=p2
                
    def operate(self, p1, p2):
        """
        return result of operation p2 to p1
        """
        assert (p1 in self.phase_lst) and (p2 in self.phase_lst)
        return self.operate_dct[(p1,p2)]
    
    def relative(self, p1, p2):
        """
        return relative position of p2 from p1
        """
        assert (p1 in self.phase_lst) and (p2 in self.phase_lst)
        return self.relative_dct[(p1,p2)]
    
    def revops(self, p):
        """
        return reverse operation of p
        """
        assert p in self.phase_lst
        return self.relative(p, "+0")
    
    def operate_int(self, p1, p2):
        assert (0<=p1<=5) and (0<=p2<=5)
        p=self.operate(self.phase_lst[p1], self.phase_lst[p2])
        return self.phase_lst.index(p)
        
    def relative_int(self, p1, p2):
        assert (0<=p1<=5) and (0<=p2<=5)
        p=self.relative(self.phase_lst[p1], self.phase_lst[p2])
        return self.phase_lst.index(p)

    def revops_int(self, p):
        assert (0 <= p <= 5)
        p = self.revops(self.phase_lst[p])
        return self.phase_lst.index(p)

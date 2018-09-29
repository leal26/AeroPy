
class structure():
    def __init__(self, geometry_parent, geometry_child, model='beam'):
        self.g_p = geometry_parent
        self.g_c = geometry_child
        self.model = model
        self.a_p = [0, 0, 0, 0]

    def displacement(self, input, x2=0, diff=None, input_type='x1'):
        parent = self.g_p.r(input, x2=x2, input_type=input_type, diff=diff)
        child = self.g_c.r(input, x2=x2, input_type=input_type, diff=diff)
        output = child - parent
        return(output)

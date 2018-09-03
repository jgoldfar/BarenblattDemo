"""Implements a logging service specialized for the pde package."""
import logging
import types


class plog(logging.getLoggerClass()):
    def __init__(self, mod=None):
        if type(mod) is types.ClassType or type(mod) is types.TypeType:
            super(plog,self).__init__('pde.'+self.get_class_text(mod))
        elif type(mod) is types.StringType:
            super(plog,self).__init__('pde.'+mod)
        else:
            raise TypeError(self.type_err('mod','class, type, or string',type(mod)))
        
        self.setLevel(logging.INFO)
        t=logging.StreamHandler()
        t.setFormatter(
                    logging.Formatter(
                                    '%(asctime)s - %(name)s in %(funcName)s - %(levelname)s\n %(message)s\n'
                                    )
                       )
        self.addHandler(t)
    
    @staticmethod    
    def bsp(num):
        """Utility method which gives a string of num backspaces"""
        return chr(8)*num
    
    @staticmethod
    def type_err(varname,tex,tgot,msg=''):
        """Standardizes the type error message"""
        return '{} of type {}, expected {}.  {}.'.format(varname,tgot,tex,msg)
    
    @staticmethod
    def implement_err(objname,fname,msg=None):
        """Standardizes the NotImplemented Error message"""
        if msg is None:
            t=''
        else:
            t=msg+': '
        
        return t+'{} does not implement {}.'.format(objname,fname)
    @staticmethod
    def newinst_info(desc, classin):
        """Standardizes the instantiation message."""
        return 'Created new {} object of type {}'.format(desc,classin)
    
    @classmethod
    def get_class_text(cls,clin,stripbase=True):
        """
        Get qualified class name information from given type.
        input:
                clin: Type or Class.
                stripbase: Boolean.  If true, remove the first entry
                            in the list making up the qualified name.
        """
        
        if type(clin) in set([types.TypeType, types.ClassType]): 
            s=str(clin)
            s=s.split("'")
            s=s[1]
            if stripbase:
                if '.' in s:
                    s=s.split('.')
                    return '.'.join(s[1:])
                        
            return s
        raise TypeError(cls.type_err('clin','class, type, or string',type(clin)))
    
"""Implements a logging service specialized for the pde package."""
import logging
import types


class plog(logging.getLoggerClass()):
    def __init__(self, mod=None):
        super(plog,self).__init__('pde.'+mod)
        
        self.setLevel(logging.WARN)
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
    
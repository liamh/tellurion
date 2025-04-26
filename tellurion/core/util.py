def ensure_1d(obj):
    '''Make object with shape `(n,)` where n=1 if obj is a scalar'''
    if obj.shape==():
        return obj.reshape(1,)
    else:
        return obj

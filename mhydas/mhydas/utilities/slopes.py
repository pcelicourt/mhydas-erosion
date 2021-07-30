from mhydas.mhydas.utilities import variablesdefinition


def vectorize(local_param):
    try:
        return [local_param[variablesdefinition.pente_1], local_param[variablesdefinition.pente_2],
                local_param[variablesdefinition.pente_3]
                ]
    except KeyError:
        raise ValueError('At least one of the slope values are not found.')

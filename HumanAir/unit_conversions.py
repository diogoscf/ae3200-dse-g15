nautical_mile = 1852 # meter
foot = 0.3048 # meter

def ft_to_m(a):
    return a * foot

def m_to_ft(a):
    return a / foot

def nm_to_m(a):
    return a * nautical_mile

def m_to_nm(a):
    return a / nautical_mile

def ft_min_to_m_s(a):
    return a * foot / 60

def m_s_to_ft_min(a):
    return a / foot * 60

def kt_to_m_s(a):
    return a * nautical_mile / 60 / 60

def m_s_to_kt(a):
    return a / nautical_mile * 60 * 60
import numpy as np

def prepare_time_arrays(t_inj, inj, t_ip, ip):    
    it0 = np.where(t_inj > 0)[0][0]
    t_inj, inj = t_inj[it0:], inj[it0:]
    i_inj_max = np.argmax(inj)
    inj_ave = np.mean(inj[20000:25000])
    inj_re = inj[i_inj_max:]
    i_inj_0 = np.where(inj_re < inj_ave)[0][0] + i_inj_max
    t_inj_0 = t_inj[i_inj_0]*1e3 + 0.05  #[ms]
    it_decay = np.argmax(ip)
    t_decay = t_ip[it_decay]*1e3-0.4
    return t_inj_0, t_decay

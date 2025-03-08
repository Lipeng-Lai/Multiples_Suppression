import numpy as np
import pylops
from pylops.optimization.basic import lsqr
from tqdm import tqdm

def adasubtraction_lsqr(data, multiples, nwin, niter=200, damp=1e-1, len1=11, len2=10):
    
    nover = (nwin[0]//2, nwin[1]//2) 
    
    nr = (nover[0]) * (data.shape[1] // nover[0])  # 计算空间方向上的样本数量
    nt = (nover[1]) * (data.shape[2] // nover[1])  # 计算时间方向上的样本数量
    
    data_sample = data[:, :nr, :nt]
    multiples_sample = multiples[:, :nr, :nt]
    
    dimsd = (nr, nt)
    nop = nwin  
    
    nwins, dims, mwin_inends, dwin_inends = pylops.signalprocessing.patch2d_design(
        dimsd, nwin, nover, nop
    )
    
    # 创建单位操作符，用于补丁化操作
    I = pylops.Identity(nwin[0]*nwin[1])  

    # 获取不使用窗函数的补丁操作符
    PatchOpH = pylops.signalprocessing.Patch2D(I, dims, dimsd, nwin, nover, nop, tapertype='none')

    # 获取使用Hanning窗函数的补丁操作符
    PatchOp  = pylops.signalprocessing.Patch2D(I, dims, dimsd, nwin, nover, nop, tapertype='hanning')
    
    nfilt = len1 # 滤波波长 原先：11, 5
    
    # 存储矩阵
    lsqr_prim_matrix = np.empty((0, nr, nt))
    lsqr_mult_matrix = np.empty((0, nr, nt))
    
    for index in tqdm(range(data.shape[0]), desc="Processing Data", ncols=100):
        data_patched = PatchOpH.H * data_sample[index].ravel()  # 对data[index]进行补丁化
        multiples_patched = PatchOpH.H * multiples_sample[index].ravel()  # 对multiples[index]进行补丁化
    
        num_patches = nwins[0] * nwins[1]  # 计算补丁的总数
        
        data_patched = np.reshape(data_patched, (num_patches, nwin[0]*nwin[1]))  # 将data_patched重新形状化，每个补丁变为一行
        multiples_patched = np.reshape(multiples_patched, (num_patches, nwin[0]*nwin[1]))  # 将multiples_patched重新形状化
        
        filters = np.zeros((num_patches, nfilt))
        primary_est = np.zeros_like(data_patched)
        multiple_est = np.zeros_like(multiples_patched)
        
        for i in range(num_patches):
            
            data_patch_i = np.reshape(data_patched[i], (nwin[0], nwin[1]))  # 将当前补丁的数据重塑为二维形状
            mult_patch_i = np.reshape(multiples_patched[i], (nwin[0], nwin[1]))  # 将当前补丁的多次波数据重塑为二维形状
            CopStack = []  # 用于存储卷积操作符的列表
            # 构建卷积操作符
            
            for j in range(nwin[0]):
                C = pylops.utils.signalprocessing.convmtx(mult_patch_i[j], nfilt)  # 对当前行应用卷积矩阵
                Cop = pylops.basicoperators.MatrixMult(C[nfilt//2:-(nfilt//2)])  # 选择卷积矩阵的中间部分，避免边缘效应
                CopStack.append(Cop)  # 将卷积操作符添加到列表中

            dataStack = data_patch_i.ravel()  # 将数据补丁展平为一维数组
            CopStack = pylops.VStack(CopStack)  # 将所有的卷积操作符堆叠成一个操作符
            
            # 通过最小二乘法求解滤波器
            # x0 = np.exp(-np.arange(nfilt)**2 / (2 * (nfilt / 5)**2))
            filt_est = lsqr(CopStack, dataStack, x0=np.zeros(nfilt), damp=damp, niter=niter)[0]  # 寻找滤波器的解
            
            # 在滤波器值过大时进行裁剪，防止在零（或接近零）的补丁中“爆炸”
            if max(abs(filt_est)) > len2: # 原先：10
                filt_est = np.ones(nfilt)  # 如果滤波器值过大，将滤波器设置为全1
            
            # 存储结果
            filters[i] = filt_est  # 存储估计的滤波器
            multiple_est[i] = (CopStack * filt_est).T  # 应用滤波器，得到估计的多次波
            primary_est[i] = (dataStack - multiple_est[i].T).T  # 从数据中减去估计的多次波，得到估计的主波
        
        # 将补丁重新拼接回去
        primary_est = PatchOp * primary_est.ravel()  # 使用补丁操作符将估计的主波数据补丁拼接回去
        multiple_est = PatchOp * multiple_est.ravel()  # 使用补丁操作符将估计的多次波数据补丁拼接回去

        # 将拼接后的数据重新塑形为补丁前的形状
        lsqr_prim = np.reshape(primary_est, (nr, nt))  # 将拼接后的主波数据重塑为原始形状
        lsqr_mult = np.reshape(multiple_est, (nr, nt))  # 将拼接后的多次波数据重塑为原始形状

        lsqr_prim_matrix = np.append(lsqr_prim_matrix, [lsqr_prim], axis=0)
        lsqr_mult_matrix = np.append(lsqr_mult_matrix, [lsqr_mult], axis=0)
        
        
        
    return lsqr_prim_matrix, lsqr_mult_matrix
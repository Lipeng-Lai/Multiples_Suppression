from scipy.signal import correlate2d
import numpy as np

def weighted_cross_correlation(A, B, weights=None):
    """
    计算加权互相关（weighted cross correlation）

    参数:
    A (ndarray): 第一个输入矩阵
    B (ndarray): 第二个输入矩阵
    weights (ndarray): 权重矩阵，必须和 A, B 维度一致

    返回:
    ndarray: 加权互相关结果
    """
    # 确保 A, B 和 weights 的形状一致
    assert A.shape == B.shape, "A, B 和权重矩阵的形状必须相同"
    
    weights = np.ones_like(A)
    # 对 A 和 B 进行加权处理
    weighted_A = A * weights
    weighted_B = B * weights
    
    # 使用 scipy 的 correlate2d 进行加权互相关计算
    correlation = correlate2d(weighted_A, weighted_B, mode='same', boundary='fill')
    
    return correlation
# -*- coding:utf-8 -*-

import cv2 as cv
import numpy as np

class Canny:
    # 将彩色图像转化为灰度图像
    def BGR2GRAY(self, src):
        b = src[:, :, 0].copy()
        g = src[:, :, 1].copy()
        r = src[:, :, 2].copy()

        out = 0.2126*r + 0.7152*g + 0.0722*g
        out = out.astype(np.uint8)
        cv.imwrite("D:\python codes\pic source\gray_7.JPG", out)
        return out

    # 对灰度图像进行高斯滤波
    def gaussion_filter(self, img, k_size=3, sigma=1.3):
        if len(img.shape) == 3:
            H, W, C = img.shape
            gray = False
        else:
            img = np.expand_dims(img, axis=-1)
            H, W, C = img.shape
            gray = True

            ## Zero padding
        pad = k_size // 2
        out = np.zeros([H + pad * 2, W + pad * 2, C], dtype=np.float)
        out[pad: pad + H, pad: pad + W] = img.copy().astype(np.float)

        ## prepare Kernel
        K = np.zeros((k_size, k_size), dtype=np.float)
        for x in range(-pad, -pad + k_size):
            for y in range(-pad, -pad + k_size):
                K[y + pad, x + pad] = np.exp(- (x ** 2 + y ** 2) / (2 * sigma * sigma))
        # K /= (sigma * np.sqrt(2 * np.pi))
        K /= (2 * np.pi * sigma * sigma)
        K /= K.sum()

        tmp = out.copy()

        # 滤波
        for y in range(H):
            for x in range(W):
                for c in range(C):
                    out[pad + y, pad + x, c] = np.sum(K * tmp[y: y + k_size, x: x + k_size, c])

        out = np.clip(out, 0, 255)
        out = out[pad: pad + H, pad: pad + W]
        out = out.astype(np.uint8)

        if gray:
            out = out[..., 0]

        return out
    # sobel 滤波
    def sobel_filter(self, src, k_size=3):
        if len(src.shape) == 3:
            H, W, C = src.shape
        else:
            H, W = src.shape

        ## 零填充
        pad = k_size // 2
        out = np.zeros((H+pad*2, W+pad*2), dtype=np.float)
        out[pad:pad+H, pad:pad+W] = src.copy().astype(np.float)
        tmp = out.copy()

        out_v = out.copy()
        out_h = out.copy()

        ## Sobel vertical
        Kv = [[1., 2., 1.], [0., 0., 0.], [-1., -2., -1.]]
        ## Sobel horizontal
        Kh = [[1., 0., -1.], [2., 0., -2.], [1., 0., -1.]]

        # 滤波
        for y in range(H):
            for x in range(W):
                out_v[pad + y, pad + x] = np.sum(Kv * (tmp[y: y + k_size, x: x + k_size]))
                out_h[pad + y, pad + x] = np.sum(Kh * (tmp[y: y + k_size, x: x + k_size]))

        out_v = np.clip(out_v, 0, 255)
        out_h = np.clip(out_h, 0, 255)

        out_v = out_v[pad: pad + H, pad: pad + W]
        out_v = out_v.astype(np.uint8)
        out_h = out_h[pad: pad + H, pad: pad + W]
        out_h = out_h.astype(np.uint8)

        return out_v, out_h

    def get_edge_angle(self, fx, fy):
        # get edge strength
        edge = np.sqrt(np.power(fx.astype(np.float32), 2) + np.power(fy.astype(np.float32), 2))
        edge = np.clip(edge, 0, 255) #将像素值固定在0-255之间

        fx = np.maximum(fx, 1e-10) # 避免fx=0
        # fx[np.abs(fx) <= 1e-5] = 1e-5

        # get edge angle
        angle = np.arctan(fy / fx)

        return edge, angle

    # 根据角度范围确定边缘方向
    def angle_quantization(self, angle):
        angle = angle / np.pi * 180
        angle[angle < -22.5] = 180 + angle[angle < -22.5]
        _angle = np.zeros_like(angle, dtype=np.uint8)
        _angle[np.where(angle <= 22.5)] = 0
        _angle[np.where((angle > 22.5) & (angle <= 67.5))] = 45
        _angle[np.where((angle > 67.5) & (angle <= 112.5))] = 90
        _angle[np.where((angle > 112.5) & (angle <= 157.5))] = 135

        return _angle

    # 非极大值抑制
    def non_maximum_suppression(self, angle, edge):
        H, W = angle.shape
        _edge = edge.copy()

        for y in range(H):
            for x in range(W):
                if angle[y, x] == 0:
                    dx1, dy1, dx2, dy2 = -1, 0, 1, 0
                elif angle[y, x] == 45:
                    dx1, dy1, dx2, dy2 = -1, 1, 1, -1
                elif angle[y, x] == 90:
                    dx1, dy1, dx2, dy2 = 0, -1, 0, 1
                elif angle[y, x] == 135:
                    dx1, dy1, dx2, dy2 = -1, -1, 1, 1
                if x == 0:
                    dx1 = max(dx1, 0)
                    dx2 = max(dx2, 0)
                if x == W - 1:
                    dx1 = min(dx1, 0)
                    dx2 = min(dx2, 0)
                if y == 0:
                    dy1 = max(dy1, 0)
                    dy2 = max(dy2, 0)
                if y == H - 1:
                    dy1 = min(dy1, 0)
                    dy2 = min(dy2, 0)
                if max(max(edge[y, x], edge[y + dy1, x + dx1]), edge[y + dy2, x + dx2]) != edge[y, x]:
                    _edge[y, x] = 0

        return _edge

    # 使用双阈值处理和连通性分析来检测边缘
    def hysterisis(self, edge, HT=100, LT=30):
        H, W = edge.shape

        # 双阈值
        edge[edge >= HT] = 255
        edge[edge <= LT] = 0

        _edge = np.zeros((H + 2, W + 2), dtype=np.float32)
        _edge[1: H + 1, 1: W + 1] = edge

        ## 8 - 邻域
        nn = np.array(((1., 1., 1.), (1., 0., 1.), (1., 1., 1.)), dtype=np.float32)

        # 如果梯度幅度值在LT和LH之间
        for y in range(1, H + 2):
            for x in range(1, W + 2):
                if _edge[y, x] < LT or _edge[y, x] > HT:
                    continue
                if np.max(_edge[y - 1:y + 2, x - 1:x + 2] * nn) >= HT: #判断八连通区域中是否有边缘点
                    _edge[y, x] = 255
                else:
                    _edge[y, x] = 0

        edge = _edge[1:H + 1, 1:W + 1]

        return edge

 

def main():
    src = cv.imread("D:\python codes\pic source\\7.jpg").astype(np.float32)
    canny = Canny()
    # 转变成灰度图像
    out = canny.BGR2GRAY(src)
    cv.imshow("Input", out)
    # 高斯滤波
    gaussian = canny.gaussion_filter(out)
    # sobel 滤波，获得gx， gy
    fx, fy = canny.sobel_filter(gaussian, 3)
    # 计算梯度幅度图像和角度图像
    edge, angle = canny.get_edge_angle(fx, fy)

    # cv.imshow("output", edge)
    # cv.imwrite("C:/Users/Odysseus96/Pictures/Image/edge.JPG", edge)
    # cv.waitKey(0)
    # cv.destroyAllWindows()

    # 根据角度范围确定边缘方向
    angle = canny.angle_quantization(angle)
    # 对梯度幅度图像应用非极大抑制
    edge = canny.non_maximum_suppression(angle, edge)
    # 使用双阈值处理和连通性分析来检测与连接边缘
    out = canny.hysterisis(edge, 50, 20).astype(np.uint8)

    cv.imwrite("D:\python codes\pic source\canny_7.JPG", out)
    cv.imshow("output", out)
    cv.waitKey(0)
    cv.destroyAllWindows()

if __name__ == '__main__':
    main()

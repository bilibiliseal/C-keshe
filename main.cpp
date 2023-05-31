#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include<string>
#include<tuple>
#include <limits>
#include <map>
using namespace std;


// 抽象命令接口
class Command {
public:
    virtual ~Command() {}
    virtual void execute() = 0;
};
// 定义点的数据结构
struct Point {
    double x;
    double y;
};
// 定义图层的数据结构
struct Layer {
    std::string name;
    std::vector<Point> points;
};
// 边的数据结构
struct Edge {
    int p1;
    int p2;
};

// 三角形的数据结构
struct Triangle {
    int p1;
    int p2;
    int p3;
};

// 接收者 - 图层管理器
class LayerManager {

public:
    void calculateConvexHull() {
        std::cout << "执行凸包计算..." << std::endl;
        // 执行凸包计算的具体实现
        std::vector<Point> convexHull = computeConvexHull(points);
        // 输出凸包点集
        std::cout << "凸包点集：" << std::endl;
        for (const Point& p : convexHull) {
            std::cout << "(" << p.x << ", " << p.y << ")" << std::endl;
        }
        
}
    

    void calculateBuffer() {
        std::cout << "执行缓冲区计算..." << std::endl;
        // 执行缓冲区计算的具体实现
         
        std::vector<Point> points = getPoints(); // 假设获取需要计算缓冲区的点集
        for (const Point& point : points) {
            // 执行缓冲区计算的操作
            double bufferedX = calculateBufferedX(point.x);
            double bufferedY = calculateBufferedY(point.y);
            std::cout << "缓冲区计算结果：(" << bufferedX << ", " << bufferedY << ")" << std::endl;
        }
    }


    void performOverlayAnalysis() {
        std::cout << "执行叠加分析..." << std::endl;
        // 两个图层 layer1 和 layer2
        Layer layer1;
        layer1.name = "Layer 1";
        layer1.points = { {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0} };

        Layer layer2;
        layer2.name = "Layer 2";
        layer2.points = { {3.0, 4.0}, {7.0, 8.0}, {9.0, 10.0} };

        // 合并两个图层的点集合并去重
        std::vector<Point> mergedPoints = layer1.points;
        mergedPoints.insert(mergedPoints.end(), layer2.points.begin(), layer2.points.end());
        std::sort(mergedPoints.begin(), mergedPoints.end(), [](const Point& p1, const Point& p2) {
            return p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y);
            });
        auto last = std::unique(mergedPoints.begin(), mergedPoints.end(), [](const Point& p1, const Point& p2) {
            return std::abs(p1.x - p2.x) < std::numeric_limits<double>::epsilon() &&
                std::abs(p1.y - p2.y) < std::numeric_limits<double>::epsilon();
            });
        mergedPoints.erase(last, mergedPoints.end());

        // 输出结果
        std::cout << "叠加分析结果：" << std::endl;
        std::cout << "图层名称: " << layer1.name << " + " << layer2.name << std::endl;
        std::cout << "合并后的点集合: " << std::endl;
        for (const Point& point : mergedPoints) {
            std::cout << "(" << point.x << ", " << point.y << ")" << std::endl;
        }
    }
    

    void calculateDelaunayTriangulation() {
        std::cout << "执行Delaunay三角网计算..." << std::endl;
        // 一个点集合
        std::vector<Point> points = { {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0} };

        // 执行Delaunay三角网计算
        std::vector<Edge> edges;
        std::vector<Triangle> triangles;
        computeDelaunayTriangulation(points, edges, triangles);

        // 输出结果
        std::cout << "Delaunay三角网结果：" << std::endl;
        std::cout << "边集合：" << std::endl;
        for (const Edge& edge : edges) {
            std::cout << edge.p1 << " -- " << edge.p2 << std::endl;
        }
        std::cout << "三角形集合：" << std::endl;
        for (const Triangle& triangle : triangles) {
            std::cout << triangle.p1 << " -- " << triangle.p2 << " -- " << triangle.p3 << std::endl;
        }
    }

    void calculateVoronoiDiagram() {
        std::cout << "执行Voronoi图计算..." << std::endl;

        // 有一组点集
        vector<Point> points = {
            {2.0, 3.0},
            {4.0, 1.0},
            {5.0, 7.0},
            {6.0, 3.0},
            {9.0, 5.0}
        };

        // 计算Voronoi图
        for (const Point& p : points) {
            // 找到与当前点p距离最近的点
            Point nearestPoint;
            double minDistance = std::numeric_limits<double>::max();

            for (const Point& q : points) {
                if (&p != &q) {  // 排除当前点p自身
                    double distance = calculateDistance(p, q);
                    if (distance < minDistance) {
                        minDistance = distance;
                        nearestPoint = q;
                    }
                }
            }

            // 输出点p的区域
            std::cout << "点 (" << p.x << ", " << p.y << ") 的区域由最近的点 (" << nearestPoint.x << ", " << nearestPoint.y << ") 组成" << std::endl;
        }
    }

    

    void calculateGeometryStatistics() {
        std::cout << "执行几何统计分析..." << std::endl;
        // 计算平均距离和标准差
        double sumDistances = 0.0;
        for (size_t i = 0; i < points.size(); ++i) {
            for (size_t j = i + 1; j < points.size(); ++j) {
                double dx = points[i].x - points[j].x;
                double dy = points[i].y - points[j].y;
                double distance = sqrt(dx * dx + dy * dy);
                sumDistances += distance;
            }
        }
        double meanDistance = sumDistances / (points.size() * (points.size() - 1) / 2);

        double sumSquaredDistances = 0.0;
        for (size_t i = 0; i < points.size(); ++i) {
            for (size_t j = i + 1; j < points.size(); ++j) {
                double dx = points[i].x - points[j].x;
                double dy = points[i].y - points[j].y;
                double distance = sqrt(dx * dx + dy * dy);
                double squaredDistance = distance * distance;
                sumSquaredDistances += squaredDistance;
            }
        }
        double variance = (sumSquaredDistances - (sumDistances * sumDistances) / (points.size() * (points.size() - 1))) / (points.size() * (points.size() - 1) / 2);
        double standardDeviation = sqrt(variance);
        cout << "平均距离: " << meanDistance << endl;
        cout << "标准差: " << standardDeviation << endl;

    }

    void calculateColorStatistics() {
        std::cout << "执行颜色统计分析..." << std::endl;
        // 统计颜色数据
        map<string, int> colorCount;
        for (const string& color : colors) {
            colorCount[color]++;
        }
        // 输出颜色统计结果
        for (const auto& pair : colorCount) {
            cout << "颜色 " << pair.first << " 的数量为 " << pair.second << endl;
        }
    }

    void addColor(const string& color) {
        colors.push_back(color);
    
    }
private:
    std::vector<Point> points = {
        {1.0, 1.0},
        {2.0, 3.0},
        {4.0, 2.0},
        {3.0, 4.0},
        {5.0, 1.0},
        {6.0, 3.0},
        {7.0, 2.0},
        {5.0, 5.0}
    };

    std::vector<Point> computeConvexHull(const std::vector<Point>& points) {
        // 具体的凸包计算算法实现
        std::vector<Point> convexHull;

        // 找到最左下的点作为起始点
        Point start = findStartPoint(points);

        // 将起始点加入凸包
        convexHull.push_back(start);

        // 根据极角排序点集
        std::vector<Point> sortedPoints = sortPointsByPolarAngle(points, start);

        // Graham扫描算法构建凸包
        for (const Point& p : sortedPoints) {
            while (convexHull.size() >= 2 && !isTurnLeft(convexHull[convexHull.size() - 2], convexHull[convexHull.size() - 1], p)) {
                convexHull.pop_back();
            }
            convexHull.push_back(p);
        }

        return convexHull;
    }

    Point findStartPoint(const std::vector<Point>& points) {
        Point start = points[0];
        for (const Point& p : points) {
            if (p.y < start.y || (p.y == start.y && p.x < start.x)) {
                start = p;
            }
        }
        return start;
    }

    std::vector<Point> sortPointsByPolarAngle(const std::vector<Point>& points, const Point& reference) {
        std::vector<Point> sortedPoints = points;
        std::sort(sortedPoints.begin(), sortedPoints.end(), [&reference](const Point& a, const Point& b) {
            double angleA = atan2(a.y - reference.y, a.x - reference.x);
            double angleB = atan2(b.y - reference.y, b.x - reference.x);
            return angleA < angleB;
            });
        return sortedPoints;
    }
    bool isTurnLeft(const Point& a, const Point& b, const Point& c) {
        double crossProduct = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        return crossProduct > 0;
    }
    // 获取需要计算缓冲区的点集的函数
    std::vector<Point> getPoints() {
        std::vector<Point> points;
        // 添加获取点集的代码

        points.push_back({ 1.0, 2.0 });
        points.push_back({ 3.0, 4.0 });
        points.push_back({ 5.0, 6.0 });
        return points;
    }

    // 计算缓冲区X坐标的函数
    double calculateBufferedX(double x) {
        // 添加计算缓冲区X坐标的代码
        
        return x + 1.0;
    }

    // 计算缓冲区Y坐标的函数
    double calculateBufferedY(double y) {
        // 添加计算缓冲区Y坐标的代码
        return y - 1.0;
    }
    void computeDelaunayTriangulation(const std::vector<Point>& points, std::vector<Edge>& edges, std::vector<Triangle>& triangles) {
        // 三角网计算的具体实现
        assert(points.size() >= 3 && "至少需要三个点才能进行三角网计算");

       //实现Delaunay三角网计算算法
      //连接了所有的点构成三角形

        for (int i = 0; i < points.size() - 2; ++i) {
            for (int j = i + 1; j < points.size() - 1; ++j) {
                for (int k = j + 1; k < points.size(); ++k) {
                    Triangle triangle;
                    triangle.p1 = i;
                    triangle.p2 = j;
                    triangle.p3 = k;
                    triangles.push_back(triangle);
                }
            }
        }

        // 构造边集合
        for (const Triangle& triangle : triangles) {
            edges.push_back({ triangle.p1, triangle.p2 });
            edges.push_back({ triangle.p2, triangle.p3 });
            edges.push_back({ triangle.p3, triangle.p1 });
        }

        // 去除重复的边
        auto isDuplicateEdge = [](const Edge& edge1, const Edge& edge2) {
            return (edge1.p1 == edge2.p1 && edge1.p2 == edge2.p2) || (edge1.p1 == edge2.p2 && edge1.p2 == edge2.p1);
        };

        std::sort(edges.begin(), edges.end(), [](const Edge& edge1, const Edge& edge2) {
            return std::tie(edge1.p1, edge1.p2) < std::tie(edge2.p1, edge2.p2);
            });
        edges.erase(std::unique(edges.begin(), edges.end(), isDuplicateEdge), edges.end());
    }
    double calculateDistance(const Point& p, const Point& q) {
        double dx = p.x - q.x;
        double dy = p.y - q.y;
        return sqrt(dx * dx + dy * dy);
    }
    vector<string> colors;  // 保存颜色数据
};

// 具体命令 - 凸包计算命令
class ConvexHullCommand : public Command {
private:
    LayerManager* layerManager;

public:
    ConvexHullCommand(LayerManager* manager) : layerManager(manager) {}

    void execute() {
        layerManager->calculateConvexHull();
    }
};

// 具体命令 - 缓冲区计算命令
class BufferCommand : public Command {
private:
    LayerManager* layerManager;

public:
    BufferCommand(LayerManager* manager) : layerManager(manager) {}

    void execute() {
        layerManager->calculateBuffer();
    }
};

// 具体命令 - 叠加分析命令
class OverlayAnalysisCommand : public Command {
private:
    LayerManager* layerManager;

public:
    OverlayAnalysisCommand(LayerManager* manager) : layerManager(manager) {}

    void execute() {
        layerManager->performOverlayAnalysis();
    }
};

// 具体命令 - Delaunay三角网计算命令
class DelaunayTriangulationCommand : public Command {
private:
    LayerManager* layerManager;

public:
    DelaunayTriangulationCommand(LayerManager* manager) : layerManager(manager) {}

    void execute() {
        layerManager->calculateDelaunayTriangulation();
    }
};

// 具体命令 - Voronoi图计算命令
class VoronoiDiagramCommand : public Command {
private:
    LayerManager* layerManager;

public:
    VoronoiDiagramCommand(LayerManager* manager) : layerManager(manager) {}
    void execute() {
        layerManager->calculateVoronoiDiagram();
    }
};

// 具体命令 - 几何统计分析命令
class GeometryStatisticsCommand : public Command {
private:
    LayerManager* layerManager;

public:
    GeometryStatisticsCommand(LayerManager* manager) : layerManager(manager) {}
    void execute() {
        layerManager->calculateGeometryStatistics();
    }
};

// 具体命令 - 颜色统计分析命令
class ColorStatisticsCommand : public Command {
private:
    LayerManager* layerManager;

public:
    ColorStatisticsCommand(LayerManager* manager) : layerManager(manager) {}
    void execute() {
        layerManager->calculateColorStatistics();
    }
};

// 请求者 - 图层分析器
class LayerAnalyzer {
private:
    std::vector<Command*> commands;

public:
    void addCommand(Command* command) {
        commands.push_back(command);
    }
    void executeCommands() {
        for (Command* command : commands) {
            command->execute();
        }
        commands.clear();
    }
};

int main() {
    // 创建图层管理器和图层分析器
    LayerManager layerManager;
    LayerAnalyzer layerAnalyzer;
    // 创建各种分析命令
    ConvexHullCommand convexHullCommand(&layerManager);
    BufferCommand bufferCommand(&layerManager);
    OverlayAnalysisCommand overlayAnalysisCommand(&layerManager);
    DelaunayTriangulationCommand delaunayTriangulationCommand(&layerManager);
    VoronoiDiagramCommand voronoiDiagramCommand(&layerManager);
    GeometryStatisticsCommand geometryStatisticsCommand(&layerManager);
    ColorStatisticsCommand colorStatisticsCommand(&layerManager);

    // 添加需要执行的命令
    layerAnalyzer.addCommand(&convexHullCommand);
    layerAnalyzer.addCommand(&bufferCommand);
    layerAnalyzer.addCommand(&overlayAnalysisCommand);
    layerAnalyzer.addCommand(&delaunayTriangulationCommand);
    layerAnalyzer.addCommand(&voronoiDiagramCommand);
    layerAnalyzer.addCommand(&geometryStatisticsCommand);
    layerAnalyzer.addCommand(&colorStatisticsCommand);
    // 添加颜色数据
    layerManager.addColor("Red");
    layerManager.addColor("Blue");
    layerManager.addColor("Green");
    layerManager.addColor("Red");
    layerManager.addColor("Yellow");
    layerManager.addColor("Blue");

    // 执行命令
    layerAnalyzer.executeCommands();

    return 0;
}


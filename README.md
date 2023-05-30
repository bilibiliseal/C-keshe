# C-keshe
课设
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
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

// 接收者 - 图层管理器
class LayerManager {
public:
    void calculateConvexHull() {
        std::cout << "执行凸包计算..." << std::endl;
        
};
    

    void calculateBuffer() {
        std::cout << "执行缓冲区计算..." << std::endl;
        // 执行缓冲区计算的具体实现
    }

    void performOverlayAnalysis() {
        std::cout << "执行叠加分析..." << std::endl;
        // 执行叠加分析的具体实现
    }

    void calculateDelaunayTriangulation() {
        std::cout << "执行Delaunay三角网计算..." << std::endl;
        // 执行Delaunay三角网计算的具体实现
    }

    void calculateVoronoiDiagram() {
        std::cout << "执行Voronoi图计算..." << std::endl;
        // 执行Voronoi图计算的具体实现
    }

    void calculateGeometryStatistics() {
        std::cout << "执行几何统计分析..." << std::endl;
        // 执行几何统计分析的具体实现
    }

    void calculateColorStatistics() {
        std::cout << "执行颜色统计分析..." << std::endl;
        // 执行颜色统计分析的具体实现
    }
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

    // 执行命令
    layerAnalyzer.executeCommands();

    return 0;
}


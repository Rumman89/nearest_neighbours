#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

struct Point3D {
    double x, y, z;
};

struct KDNode {
    Point3D point;
    KDNode* left;
    KDNode* right;
};

class KDTree {
public:
    KDTree(const std::vector<Point3D>& points) {
        root = buildTree(points, 0);
    }

    double findClosestDistance(const Point3D& target) {
        return findClosestDistance(root, target, 0);
    }

private:
    KDNode* root;

    KDNode* buildTree(std::vector<Point3D> points, int depth) {
        if (points.empty()) return nullptr;

        int axis = depth % 3;
        std::sort(points.begin(), points.end(), [axis](const Point3D& a, const Point3D& b) {
            return axis == 0 ? a.x < b.x : axis == 1 ? a.y < b.y : a.z < b.z;
        });

        int medianIndex = points.size() / 2;
        KDNode* node = new KDNode{ points[medianIndex], nullptr, nullptr };
        std::vector<Point3D> leftPoints(points.begin(), points.begin() + medianIndex);
        std::vector<Point3D> rightPoints(points.begin() + medianIndex + 1, points.end());

        node->left = buildTree(leftPoints, depth + 1);
        node->right = buildTree(rightPoints, depth + 1);

        return node;
    }

    double findClosestDistance(KDNode* node, const Point3D& target, int depth) {
        if (!node) return std::numeric_limits<double>::max();

        double distance = euclideanDistance(node->point, target);
        int axis = depth % 3;

        KDNode* nextNode = nullptr;
        KDNode* otherNode = nullptr;
        if (axis == 0) {
            nextNode = target.x < node->point.x ? node->left : node->right;
            otherNode = target.x < node->point.x ? node->right : node->left;
        }
        else if (axis == 1) {
            nextNode = target.y < node->point.y ? node->left : node->right;
            otherNode = target.y < node->point.y ? node->right : node->left;
        }
        else {
            nextNode = target.z < node->point.z ? node->left : node->right;
            otherNode = target.z < node->point.z ? node->right : node->left;
        }

        double closestDistance = findClosestDistance(nextNode, target, depth + 1);
        if (distance < closestDistance) closestDistance = distance;

        double axisDistance = axis == 0 ? std::abs(target.x - node->point.x) :
            axis == 1 ? std::abs(target.y - node->point.y) :
            std::abs(target.z - node->point.z);
        if (axisDistance < closestDistance) {
            double otherDistance = findClosestDistance(otherNode, target, depth + 1);
            if (otherDistance < closestDistance) closestDistance = otherDistance;
        }

        return closestDistance;
    }

    double euclideanDistance(const Point3D& a, const Point3D& b) {
        return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2) + std::pow(a.z - b.z, 2));
    }
};

class CoverTree {
public:
    CoverTree(const std::vector<Point3D>& points) {
        root = new Node(points[0]);
        for (size_t i = 1; i < points.size(); ++i) {
            insert(root, points[i]);
        }
    }

    double findClosestDistance(const Point3D& target) {
        return findClosestDistance(root, target);
    }

private:
    struct Node {
        Point3D point;
        std::vector<Node*> children;

        Node(const Point3D& p) : point(p) {}
    };

    Node* root;

    void insert(Node* node, const Point3D& point) {
        double distance = euclideanDistance(node->point, point);
        for (Node* child : node->children) {
            if (distance <= euclideanDistance(child->point, point)) {
                insert(child, point);
                return;
            }
        }
        node->children.push_back(new Node(point));
    }

    double findClosestDistance(Node* node, const Point3D& target) {
        double closestDistance = euclideanDistance(node->point, target);
        for (Node* child : node->children) {
            double childDistance = findClosestDistance(child, target);
            if (childDistance < closestDistance) {
                closestDistance = childDistance;
            }
        }
        return closestDistance;
    }

    double euclideanDistance(const Point3D& a, const Point3D& b) {
        return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2) + std::pow(a.z - b.z, 2));
    }
};


int main() {
    std::vector<Point3D> points = { {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0} };
    KDTree kd_tree(points);

    Point3D target = { 2.0, 3.0, 4.0 };
    double closestDistance = kd_tree.findClosestDistance(target);

    std::cout << "Closest distance: " << closestDistance << std::endl;

    CoverTree cover_tree(points);
    closestDistance = cover_tree.findClosestDistance(target);

    std::cout << "Closest distance: " << closestDistance << std::endl;

    return 0;
}

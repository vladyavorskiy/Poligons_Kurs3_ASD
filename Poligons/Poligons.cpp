#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm> 
#include <random>

using namespace  std;


class Point {
public:
    double x, y;

    Point(double x, double y) : x(x), y(y){}


};

class polygon {
private:
    vector<shared_ptr<Point>> vertices;
    bool isConvex;

public:
    // Конструктор
    polygon(bool convex = true) : isConvex(convex) {}

    // Деструктор
    ~polygon() = default;

    vector<pair<shared_ptr<Point>, shared_ptr<Point>>> findPossibleInsertionPoints(double newX, double newY) {
        vector<pair<shared_ptr<Point>, shared_ptr<Point>>> possibleInsertionPoints;
        int count = vertices.size();
        for (int i = 0; i < count; i++) {
            int j = (i + 1) % count;
            bool flag = true;
            for (int k1 = 0; k1 < count; k1++) {
                int l1 = (k1 + 1) % count;
                if ((k1 == i) || (l1 == i)) {
                    continue;
                }
                if (isIntersect(newX, newY, vertices[i]->x, vertices[i]->y, vertices[k1]->x, vertices[k1]->y, vertices[l1]->x, vertices[l1]->y)) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                for (int k2 = 0; k2 < count; k2++) {
                    int l2 = (k2 + 1) % count;
                    if ((k2 == j) || (l2 == j)) {
                        continue;
                    }
                    if (isIntersect(newX, newY, vertices[j]->x, vertices[j]->y, vertices[k2]->x, vertices[k2]->y, vertices[l2]->x, vertices[l2]->y)) {
                        flag = false;
                        break;
                    }
                }
            }
            if (flag) {
                possibleInsertionPoints.emplace_back(vertices[i], vertices[j]);
            }
        }
        return possibleInsertionPoints;
    }


    bool isIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
        //x1=newx1 y1=newy1 //x2=i.x  y2=i.y
        //x3=k.x  y3=k.y //x4=l.x  y4=l.y
        double verx1 = ((x3 - x1) * (y2 - y1)) - ((y3 - y1) * (x2 -x1));
        double verx2 = ((x3 - x1) * (y4 - y3)) - ((y3 - y1) * (x4 - x3));
        double niz = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1));

        if ((verx1 / niz) >= 0 and (verx1 / niz) <= 1) {
            if ((verx2 / niz) >= 0 and (verx2 / niz) <= 1) {
                return true;
            }
        }
        return false;
    }


    void addVertex(double x, double y) {
        int count = vertices.size();
        if (count == 0) {
            vertices.push_back(make_shared<Point>(x, y));
        }
        else if (count < 3) {
            vertices.push_back(make_shared<Point>(x, y));
        }
        else {
            vector<pair<shared_ptr<Point>, shared_ptr<Point>>> possiblePoints = findPossibleInsertionPoints(x, y);
            int i;
            if (possiblePoints.size() == 0) {
                cout << "Невозможно вставить точку с координатами (" << x << ", " << y << ") : "<< endl;
                return;
            }
            else if (possiblePoints.size() == 1) {
                i = 0;
            }
            else {
                cout << "Возможное положение новой точки (" << x << ", " << y << "):" << endl;
                for (int i = 0; i < possiblePoints.size(); ++i) {
                    cout << i << ". " << "Между вершинами: " << "(" << possiblePoints[i].first->x << ", " << possiblePoints[i].first->y << ") и (" << possiblePoints[i].second->x << ", " << possiblePoints[i].second->y << ")" << endl;
                }
                cin >> i;
            }
            

            auto insertPosition = find(vertices.begin(), vertices.end(), possiblePoints[i].first);
            int indexPos = distance(vertices.begin(), insertPosition);
            vertices.insert(insertPosition + 1, make_shared<Point>(x, y));
            
        }
        return;
    }


    void addVertexRandom(double x, double y) {
        int count = vertices.size();
        if (count == 0) {
            vertices.push_back(make_shared<Point>(x, y));
        }
        else if (count < 3) {
            vertices.push_back(make_shared<Point>(x, y));
            auto firstIndex = 0;
            auto lastIndex = count - 1;
        }
        else {
            vector<pair<shared_ptr<Point>, shared_ptr<Point>>> possiblePoints = findPossibleInsertionPoints(x, y);

            random_device rd;
            mt19937 generator(rd());

            uniform_int_distribution<int> distribution(0, possiblePoints.size() - 1);
            int i = distribution(generator);

            auto insertPosition = find(vertices.begin(), vertices.end(), possiblePoints[i].first);
            int indexPos = distance(vertices.begin(), insertPosition);
            vertices.insert(insertPosition + 1, make_shared<Point>(x, y));



        }
        return;
    }

    // Удаление вершины по индексу
    void removeVertex(int index) {
        int count = vertices.size();
        bool flag = true;
            if (index >= 0 && index < vertices.size()) {
                int i = (index + count - 1) % count;
                int j = (index + 1) % count;
                for (int k = 0; k < count; k++) {
                    int l = (k + 1) % count;
                    if ((k == i) || (l == i) || (k == j) || (l == j)) {
                        continue;
                    }
                    if (isIntersect(vertices[i]->x, vertices[i]->y, vertices[j]->x, vertices[j]->y, vertices[k]->x, vertices[k]->y, vertices[l]->x, vertices[l]->y)) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    vertices.erase(vertices.begin() + index);
                    return;
                }
                else {
                    cout << "Невозможно удалить вершину с таким индексом" << endl;
                    return;
                }
            }
            else {
                cout << "Вершины с таким индексом нет" << endl;
            }
    }

    // Удаление вершины по координатам
    void removeVertex(double x, double y) {
        int count = vertices.size();
        bool flag = true;
        for (int n = 0; n < count; n++) {
            if (vertices[n]->x == x && vertices[n]->y == y) {
                int i = (n + count - 1) % count;
                int j = (n + 1) % count;
                for (int k = 0; k < count; k++) {
                    int l = (k + 1) % count;
                    if ((k == i) || (l == i) || (k == j) || (l == j)) {
                        continue;
                    }
                    if (isIntersect(vertices[i]->x, vertices[i]->y, vertices[j]->x, vertices[j]->y, vertices[k]->x, vertices[k]->y, vertices[l]->x, vertices[l]->y)) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    vertices.erase(vertices.begin() + n);
                    return;
                }
                else {
                    cout << "Невозможно удалить вершину с такими координатами" << endl;
                    return;
                }
            }
        }
        cout << "Вершины с такими координатами нет" << endl;
    }

    // Метод для нахождения площади многоугольника
    double calculateArea() const {
        if (!isConvex) {
            std::cerr << "Cannot calculate area for a non-convex polygon." << std::endl;
            return 0.0;
        }

        double area = 0.0;
        int n = vertices.size();

        for (int i = 0; i < n; ++i) {
            int j = (i + 1) % n;
            area += (vertices[i]->x * vertices[j]->y - vertices[j]->x * vertices[i]->y);
        }

        /*for (int i = 0; i < n; ++i) {
            int next = (i + 1) % n;
            int prev = (i - 1) % n
            area += (vertices[i]->x * vertices[j]->y - vertices[j]->x * vertices[i]->y);
        }*/

        return abs(area) / 2.0;
    }

    // Метод для нахождения периметра многоугольника
    double calculatePerimeter() const {
        double perimeter = 0.0;
        int n = vertices.size();

        for (int i = 0; i < n; ++i) {
            int j = (i + 1) % n;
            double dx = vertices[j]->x - vertices[i]->x;
            double dy = vertices[j]->y - vertices[i]->y;
            perimeter += sqrt(dx * dx + dy * dy);
        }

        return perimeter;
    }

    // Метод для вывода многоугольника на экран
    void printPolygonVertex() const {
        cout << "Polygon vertices:" << endl;
        int count = vertices.size();
        if (count == 0) {
            cout << "Polygon no:" << endl;
            return;
        }

        for (int i = 0; i < count;i++) {
            cout << i << ": (" << vertices[i]->x << ", " << vertices[i]->y << ")" << endl;
        }

        cout << "For Desmos:" << endl;
        cout << "polygon(";
        for (int i = 0; i < count - 1; ++i) {
            cout << "(" << vertices[i]->x << ", " << vertices[i]->y << "), ";
        }
        cout << "(" << vertices[count-1]->x << ", " << vertices[count-1]->y << "))" << endl << endl;
    }

   
    bool pointInside(double x, double y) {
        int count = vertices.size();
        bool inside = false;
        for (int i = 0; i < count; i++) {
            int j = (i + 1) % count;
            if (x == vertices[i]->x && y == vertices[i]->y) {
                return true;
            }
            if (((vertices[i]->y > y) != (vertices[j]->y > y)) && 
                (x < vertices[i]->x + (vertices[j]->x - vertices[i]->x) * (y - vertices[i]->y) / (vertices[j]->y - vertices[i]->y)))
                inside = not inside;
        }
        return inside;
    }


    vector<shared_ptr<Point>> pointIntersect(polygon polySec) {
        vector<shared_ptr<Point>> intersectVertices;

        int count_1 = vertices.size();
        int count_2 = polySec.vertices.size();
        for (int i1 = 0; i1 < count_1; i1++) {
            int j2 = (i1 + 1) % count_1;
            for (int i3 = 0; i3 < count_2; i3++) {
                int j4 = (i3 + 1) % count_2;
                if (isIntersect(vertices[i1]->x, vertices[i1]->y, vertices[j2]->x, vertices[j2]->y, polySec.vertices[i3]->x, polySec.vertices[i3]->y, polySec.vertices[j4]->x, polySec.vertices[j4]->y)) {
                    if (vertices[i1]->x == polySec.vertices[i3]->x && vertices[i1]->y == polySec.vertices[i3]->y) {
                        intersectVertices.push_back(make_shared<Point>(vertices[i1]->x, vertices[i1]->y));
                    }
                    else if (vertices[i1]->x == polySec.vertices[j4]->x && vertices[i1]->y == polySec.vertices[j4]->y) {
                        intersectVertices.push_back(make_shared<Point>(vertices[i1]->x, vertices[i1]->y));
                    }
                    else if (vertices[j2]->x == polySec.vertices[j4]->x && vertices[j2]->y == polySec.vertices[j4]->y) {
                        intersectVertices.push_back(make_shared<Point>(vertices[j2]->x, vertices[j2]->y));
                    }
                    else if (vertices[j2]->x == polySec.vertices[i3]->x && vertices[j2]->y == polySec.vertices[i3]->y) {
                        intersectVertices.push_back(make_shared<Point>(vertices[j2]->x, vertices[j2]->y));
                    }
                    else {
                        double m1 = (vertices[j2]->y - vertices[i1]->y) / (vertices[j2]->x - vertices[i1]->x);
                        double m2 = (polySec.vertices[j4]->y - polySec.vertices[i3]->y) / (polySec.vertices[j4]->x - polySec.vertices[i3]->x);

                        double b1;
                        double b2;
                        double x;
                        double y;
                        if (isinf(m1)) {
                            b2 = polySec.vertices[i3]->y - m2 * polySec.vertices[i3]->x;
                            x = vertices[i1]->x;
                            y = x * m2 + b2;
                        }
                        else if (isinf(m2)) {
                            b1 = vertices[i1]->y - m1 * vertices[i1]->x;
                            x = polySec.vertices[i3]->x;
                            y = x * m1 + b1;
                        }
                        else {
                            b1 = vertices[i1]->y - m1 * vertices[i1]->x;
                            b2 = polySec.vertices[i3]->y - m2 * polySec.vertices[i3]->x;

                            x = (b1 - b2) / (m2 - m1);
                            y = x * m1 + b1;
                        }

                        intersectVertices.push_back(make_shared<Point>(x, y));
                        //cout << "(" << x << "; " << y << ")" << endl;
                    }
                }

            }
        }
        return intersectVertices;
    }


    polygon intersectionPolygons(polygon polySec) {
        polygon intersectionPolygon(true);
        int count_1 = vertices.size();
        int count_2 = polySec.vertices.size();

        for (int i = 0; i < count_1; i++) {
            if (polySec.pointInside(vertices[i]->x, vertices[i]->y)) {
                intersectionPolygon.vertices.push_back(vertices[i]);
            }
        }
        for (int j = 0; j < count_2; j++) {
            if (pointInside(polySec.vertices[j]->x, polySec.vertices[j]->y)) {
                intersectionPolygon.vertices.push_back(polySec.vertices[j]);
            }
        }

        if (intersectionPolygon.vertices.size() == 0) {
            return intersectionPolygon;
        }
        vector<shared_ptr<Point>> intersectVertices = pointIntersect(polySec);
        
        intersectionPolygon.vertices.insert(intersectionPolygon.vertices.end(), intersectVertices.begin(), intersectVertices.end());
        
        intersectionPolygon.removeDuplicates();

        intersectionPolygon.orderPointsClockwise();

        return intersectionPolygon;
    }

   
    void removeDuplicates() {
        // Сортируем вектор, чтобы соседние дубликаты стали смежными
        sort(vertices.begin(), vertices.end(),
            [](const shared_ptr<Point>& a, const shared_ptr<Point>& b) {
                return a->x < b->x || (a->x == b->x && a->y < b->y);
            });

        // Используем std::unique для перемещения дубликатов в конец вектора
        auto last = unique(vertices.begin(), vertices.end(),
            [](const shared_ptr<Point>& a, const shared_ptr<Point>& b) {
                return a->x == b->x && a->y == b->y;
            });

        // Удаляем дубликаты из вектора
        vertices.erase(last, vertices.end());
    }


    polygon unionPolygons(polygon polySec) {
        polygon unionPolygon(true);
        int count_1 = vertices.size();
        int count_2 = polySec.vertices.size();

        for (int i = 0; i < count_1; i++) {
            if (polySec.pointInside(vertices[i]->x, vertices[i]->y) == false) {
                unionPolygon.vertices.push_back(vertices[i]);
            }
        }
        for (int j = 0; j < count_2; j++) {
            if (pointInside(polySec.vertices[j]->x, polySec.vertices[j]->y) == false) {
                unionPolygon.vertices.push_back(polySec.vertices[j]);
            }
        }

        if (unionPolygon.vertices.size() == count_1 + count_2) {
            unionPolygon.vertices.clear();
            return unionPolygon;
        }
        vector<shared_ptr<Point>> intersectVertices = pointIntersect(polySec);
        if (intersectVertices.size() == 2 &&
            intersectVertices[0]->x == intersectVertices[1]->x &&
            intersectVertices[0]->y == intersectVertices[1]->y) {
            unionPolygon.vertices.clear();
            return unionPolygon;
        }

        unionPolygon.vertices.insert(unionPolygon.vertices.end(), intersectVertices.begin(), intersectVertices.end());

        unionPolygon.removeDuplicates();



        unionPolygon.orderPointsClockwise();

        return unionPolygon;
    }


    double calculateAngle(const shared_ptr<Point>& point, const shared_ptr<Point>& center) {
        return atan2(point->y - center->y, point->x - center->x);
    }

    shared_ptr<Point> findPolygonCenter() {
        double centerX = 0.0, centerY = 0.0;
        int n = vertices.size();

        for (auto& vertex : vertices) {
            centerX += vertex->x;
            centerY += vertex->y;
        }

        centerX /= n;
        centerY /= n;

        return make_shared<Point>(centerX, centerY);
    }

    void orderPointsClockwise() {
        shared_ptr<Point> center = findPolygonCenter();

        std::sort(vertices.begin(), vertices.end(),
            [this, center](const auto& a, const auto& b) {
                return calculateAngle(a, center) < calculateAngle(b, center);
            });
    }



};





int main() {
    // Пример использования класса
    system("chcp 1251 > null");

    //polygon poly(true);

    /*poly.addVertexRandom(0.0, 0.0);
    poly.addVertexRandom(0.0, 5.0);
    poly.addVertexRandom(5.0, 5.0);
    poly.addVertexRandom(7, 6);

    poly.printPolygonVertex();
 
    poly.addVertexRandom(4, 2);

    poly.printPolygonVertex();

    poly.removeVertex(0,5);

    poly.printPolygonVertex();*/


    /*poly.addVertex(0.0, 0.0);
    poly.addVertex(0.0, 5.0);
    poly.addVertex(5.0, 5.0);
    poly.addVertex(7, 6);

    poly.printPolygonVertex();

    poly.addVertex(4, 2);

    poly.printPolygonVertex();

    poly.removeVertex(0,5);

    poly.printPolygonVertex();*/

    /*polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(1, 5);
    poly.addVertex(6, 6);
    poly.addVertex(7, 2);
    poly.addVertex(5, 0);

    poly.printPolygonVertex();

    polygon poly2(true);
    poly2.addVertex(3,2);
    poly2.addVertex(7,8);
    poly2.addVertex(12,5);
    poly2.addVertex(10,0);

    poly2.printPolygonVertex();


    poly.pointIntersect(poly2);*/



    /*polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(0, 2);
    poly.addVertex(2, 2);
    poly.addVertex(2, 0);

    poly.printPolygonVertex();

    polygon poly2(true);
    poly2.addVertex(1, 1);
    poly2.addVertex(0.5, 3);
    poly2.addVertex(3, 3);
    poly2.addVertex(3, 1);

    poly2.printPolygonVertex();*/


    //poly.pointIntersect(poly2);



    /*polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(0, 4);
    poly.addVertex(4, 4);
    poly.addVertex(4, 0);

    poly.printPolygonVertex();

    polygon poly2(true);
    poly2.addVertex(1, 2);
    poly2.addVertex(3, 6);
    poly2.addVertex(7, 6);
    poly2.addVertex(7, 0);

    poly2.printPolygonVertex();*/


    

    /*polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(0, 4);
    poly.addVertex(4, 4);
    poly.addVertex(4, 0);

    poly.printPolygonVertex();

    polygon poly2(true);
    poly2.addVertex(1, 2);
    poly2.addVertex(2, 1);
    poly2.addVertex(2, 2);
    poly2.addVertex(1, 1);

    poly2.printPolygonVertex();*/

    polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(0, 4);
    poly.addVertex(4, 4);
    poly.addVertex(4, 0);

    poly.printPolygonVertex();

    polygon poly2(true);
    poly2.addVertex(4, 4);
    poly2.addVertex(4, 8);
    poly2.addVertex(8, 8);
    poly2.addVertex(8, 4);

    poly2.printPolygonVertex();


    /*polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(0, 4);
    poly.addVertex(4, 4);
    poly.addVertex(4, 0);

    poly.printPolygonVertex();
    cout << poly.findPolygonCenter().x << " " << poly.findPolygonCenter().y << endl;

    polygon poly2(true);
    poly2.addVertex(2, 2);
    poly2.addVertex(4, 4);
    poly2.addVertex(6,4);
    poly2.addVertex(6, 2);

    poly2.printPolygonVertex();
    cout << poly2.findPolygonCenter().x << " " << poly2.findPolygonCenter().y << endl;


    polygon poly3(true);
    poly3.addVertex(1, 2);
    poly3.addVertex(3, 6);
    poly3.addVertex(7, 6);
    poly3.addVertex(7, 0);
    poly3.addVertex(3, 3);

    poly3.printPolygonVertex();
    cout << poly3.findPolygonCenter().x << " " << poly3.findPolygonCenter().y << endl;*/

    
    /*polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(0, 4);
    poly.addVertex(4, 4);
    poly.addVertex(4, 0);

    poly.printPolygonVertex();*/

    //polygon poly2(true);
    //poly2.addVertex(1, 2);
    //poly2.addVertex(0.5, 6);
    ////poly2.addVertex(5, 3);
    //poly2.addVertex(7, 0);

    /*polygon poly2(true);
    poly2.addVertex(2, 2);
    poly2.addVertex(4, 4);
    poly2.addVertex(6, 2);
    poly2.addVertex(3, 1);*/

    /*polygon poly2(true);
    poly2.addVertex(4, 0);
    poly2.addVertex(4, 4);
    poly2.addVertex(6, 4);
    poly2.addVertex(6, 0);

    poly2.printPolygonVertex();*/

    //polygon poly(true);
    //poly.addVertex(0, 0);
    //poly.addVertex(0, 6);
    //poly.addVertex(3, 6);
    ////poly.addVertex(5, 3);
    //poly.addVertex(3, 3);
    //poly.addVertex(4, 0);

    //cout << poly.calculateArea() << endl;
    //cout << poly.calculatePerimeter() << endl;
    //cout << endl;
    //
    //poly.printPolygonVertex();
    //cout << endl;

    //polygon poly2(true);
    //poly2.addVertex(1.2, 3);
    //poly2.addVertex(3, 6);
    //poly2.addVertex(7, 6);
    //poly2.addVertex(7, 0);
    //poly2.addVertex(4, 0);

    //cout << poly2.calculateArea() << endl;
    //cout << poly2.calculatePerimeter() << endl;
    //cout << endl;

    //poly2.printPolygonVertex();
    //cout << endl;


    /*polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(0, 4);
    poly.addVertex(4, 4);
    poly.addVertex(4, 0);

    cout << poly.calculateArea() << endl;
    cout << poly.calculatePerimeter() << endl;
    cout << endl;

    poly.printPolygonVertex();
    cout << endl;

    polygon poly2(true);
    poly2.addVertex(6, 0);
    poly2.addVertex(6, 4);
    poly2.addVertex(10, 4);
    poly2.addVertex(10, 0);

    cout << poly2.calculateArea() << endl;
    cout << poly2.calculatePerimeter() << endl;
    cout << endl;

    poly2.printPolygonVertex();
    cout << endl;*/



    /*polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(0, 4);
    poly.addVertex(4, 4);
    poly.addVertex(4, 0);

    poly.printPolygonVertex();

    polygon poly2(true);
    poly2.addVertex(1, 2);
    poly2.addVertex(3, 6);
    poly2.addVertex(7, 6);
    poly2.addVertex(7, 0);
    poly2.addVertex(3, 3);

    poly2.printPolygonVertex();*/



    polygon polyInter(true);

    polyInter = poly.intersectionPolygons(poly2);
    
    //cout << polyInter.findPolygonCenter()->x << " " << polyInter.findPolygonCenter()->y << endl;
    
    cout << polyInter.calculateArea() << endl;
    cout << polyInter.calculatePerimeter() << endl;
    cout << endl;
    polyInter.printPolygonVertex();
    cout << "max" << polyInter.findPolygonCenter()->x << " " << polyInter.findPolygonCenter()->y << endl;
    cout << endl;

    polygon polyUni(true);

    polyUni = poly.unionPolygons(poly2);
    cout << polyUni.findPolygonCenter()->x << " " << polyUni.findPolygonCenter()->y << endl;

    cout << polyUni.calculateArea() << endl;
    cout << polyUni.calculatePerimeter() << endl;
    cout << endl;

    polyUni.printPolygonVertex();
    cout << endl;






    //polygon poly(true);
    //poly.addVertex(0, 0);
    //poly.addVertex(0, 6);
    //poly.addVertex(3, 6);
    ////poly.addVertex(5, 3);
    //poly.addVertex(3, 3);
    //poly.addVertex(4, 0);
    //poly.addVertex(7, 8);
    //
    //poly.printPolygonVertex();

    //poly.removeVertex(5);


    //poly.printPolygonVertex();



    /*double x1, y1, x2, y2, x3, y3, x4, y4;
    polygon poly(true);

    x1 = 2;
    y1 = 2;
    x2 = 4;
    y2 = 4;

    x3 = 4;
    y3 = 4;
    x4 = 4;
    y4 = 0;

    cout << poly.isIntersect(x1, y1, x2, y2, x3, y3, x4, y4) << endl;
    vector<shared_ptr<Point>> ppp = poly.onePointIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
    cout << ppp[0]->x << " " << ppp[0]->y << endl;*/



    /*cout << "Area: " << poly.calculateArea() << std::endl;
    cout << "Perimeter: " << poly.calculatePerimeter() << std::endl;*/



    /*cout << poly.pointInside(1, 1) << endl;
    cout << poly.pointInside(2, 2) << endl;
    cout << poly.pointInside(2, 0) << endl;
    cout << poly.pointInside(-1, 1) << endl;
    cout << poly.pointInside(3, 1) << endl;*/


    return 0;
}
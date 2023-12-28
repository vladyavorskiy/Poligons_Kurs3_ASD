#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm> 

#include <random>

using namespace  std;


class Point {
public:
    double x, y;
    shared_ptr<Point> next;
    shared_ptr<Point> prev;
    Point(double x, double y) : x(x), y(y), next(nullptr), prev(nullptr){}

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
            int j;
            if (i == count - 1) {
                j = 0;
            }
            else {
                j = i + 1;
            }

            if (isNotIntersect(i, newX, newY)) {
                if (isNotIntersect(j, newX, newY)) {
                    possibleInsertionPoints.emplace_back(vertices[i], vertices[j]);
                }
            }
        }
        return possibleInsertionPoints;
    }

    bool isNotIntersect(int i, double newx1, double newy1) {
        int count = vertices.size();
        for (int k = 0; k < count; k++) {
            int l;
            if (k == count - 1) {
                l = 0;
            }
            else {
                l = k + 1;
            }
            if ((k == i) || (l == i)) {
                continue;
            }
            //x1=newx1 y1=newy1 //x2=i.x  y2=i.y
            //x3=k.x  y3=k.y //x4=l.x  y4=l.y
            double verx1 = ((vertices[k]->x - newx1) * (vertices[i]->y - newy1)) - ((vertices[k]->y - newy1) * (vertices[i]->x - newx1));
            double verx2 = ((vertices[k]->x - newx1) * (vertices[l]->y - vertices[k]->y)) - ((vertices[k]->y - newy1) * (vertices[l]->x - vertices[k]->x));
            double niz = ((vertices[l]->y - vertices[k]->y) * (vertices[i]->x - newx1)) - ((vertices[l]->x - vertices[k]->x) * (vertices[i]->y - newy1));

            if ((verx1 / niz) >= 0 and (verx1 / niz) <= 1) {
                if ((verx2 / niz) >= 0 and (verx2 / niz) <= 1) {
                    return false;
                }
            }
        }
        return true;
    }


    bool isNotIntersect(int i, int j) {
        int count = vertices.size();
        for (int k = 0; k < count; k++) {
            int l;
            if (k == count - 1) {
                l = 0;
            }
            else {
                l = k + 1;
            }
            if ((k == i) || (l == i) || (k == j) || (l == j)) {
                continue;
            }
            //x1=newx1 y1=newy1 //x2=i.x  y2=i.y
            //x3=k.x  y3=k.y //x4=l.x  y4=l.y
            double verx1 = ((vertices[k]->x - vertices[j]->x) * (vertices[i]->y - vertices[j]->y)) - ((vertices[k]->y - vertices[j]->y) * (vertices[i]->x - vertices[j]->x));
            double verx2 = ((vertices[k]->x - vertices[j]->x) * (vertices[l]->y - vertices[k]->y)) - ((vertices[k]->y - vertices[j]->y) * (vertices[l]->x - vertices[k]->x));
            double niz = ((vertices[l]->y - vertices[k]->y) * (vertices[i]->x - vertices[j]->x)) - ((vertices[l]->x - vertices[k]->x) * (vertices[i]->y - vertices[j]->y));

            if ((verx1 / niz) >= 0 and (verx1 / niz) <= 1) {
                if ((verx2 / niz) >= 0 and (verx2 / niz) <= 1) {
                    return false;
                }
            }
        }
        return true;
    }


    bool isIntersect(int i1, int j2, int x3, int y3, int x4, int y4) {
        //x1=newx1 y1=newy1 //x2=i.x  y2=i.y
        //x3=k.x  y3=k.y //x4=l.x  y4=l.y
        double verx1 = ((x3 - vertices[i1]->x) * (vertices[j2]->y - vertices[i1]->y)) - ((y3 - vertices[i1]->y) * (vertices[j2]->x - vertices[i1]->x));
        double verx2 = ((x3 - vertices[i1]->x) * (y4 - y3)) - ((y3 - vertices[i1]->y) * (x4 - x3));
        double niz = ((y4 - y3) * (vertices[j2]->x - vertices[i1]->x)) - ((x4 - x3) * (vertices[j2]->y - vertices[i1]->y));

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
            auto firstIndex = 0;
            auto lastIndex = count - 1;
            vertices[count]->next = vertices[firstIndex];
            vertices[count]->prev = vertices[lastIndex];
            vertices[lastIndex]->next = vertices[count];
            vertices[firstIndex]->prev = vertices[count];
        }
        else {
            vector<pair<shared_ptr<Point>, shared_ptr<Point>>> possiblePoints = findPossibleInsertionPoints(x, y);
            int i;
            if (possiblePoints.size() == 0) {
                cout << "imposible" << endl;
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

            vertices[(indexPos + 1) % (count + 1)]->next = vertices[(indexPos + 2) % (count + 1)];
            vertices[(indexPos + 1) % (count + 1)]->prev = vertices[indexPos];
            vertices[(indexPos)]->next = vertices[(indexPos + 1) % (count + 1)];
            vertices[(indexPos + 2) % (count + 1)]->prev = vertices[(indexPos + 1) % (count + 1)];
            
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
            vertices[count]->next = vertices[firstIndex];
            vertices[count]->prev = vertices[lastIndex];
            vertices[lastIndex]->next = vertices[count];
            vertices[firstIndex]->prev = vertices[count];
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

            vertices[(indexPos + 1) % (count + 1)]->next = vertices[(indexPos + 2) % (count + 1)];
            vertices[(indexPos + 1) % (count + 1)]->prev = vertices[indexPos];
            vertices[(indexPos)]->next = vertices[(indexPos + 1) % (count + 1)];
            vertices[(indexPos + 2) % (count + 1)]->prev = vertices[(indexPos + 1) % (count + 1)];

        }
        return;
    }

    // Удаление вершины по индексу
    void removeVertex(int index) {
        int count = vertices.size();
        if (index >= 0 && index < vertices.size()) {
            if (isNotIntersect((index - 1) % count, (index + 1) % count)) {
                vertices[(index + 1) % count]->prev = vertices[index % count]->prev;
                vertices[(index - 1) % count]->next = vertices[index % count]->next;

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
        for (int i = 0; i < count; i++) {
            if (vertices[i]->x == x && vertices[i]->y == y) {
                if (isNotIntersect((i - 1) % count, (i + 1) % count)) {
                    vertices[(i + 1) % count]->prev = vertices[i % count]->prev;
                    vertices[(i - 1) % count]->next = vertices[i % count]->next;

                    vertices.erase(vertices.begin() + i);
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
        for (auto vertex : vertices) {
            cout << "(" << vertex->x << ", " << vertex->y << ")" << endl;
        }
        int count = vertices.size();
        cout << "For Desmos:" << endl;
        cout << "polygon(";
        for (int i = 0; i < count - 1; ++i) {
            cout << "(" << vertices[i]->x << ", " << vertices[i]->y << "), ";
        }
        cout << "(" << vertices[count-1]->x << ", " << vertices[count-1]->y << "))" << endl;
    }

   
    bool pointInside(double x, double y) {
        int count = vertices.size();
        bool inside = false;
        for (int i = 0; i < count; i++) {
            int j = (i + 1) % count;
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
                if (isIntersect(i1, j2, polySec.vertices[i3]->x, polySec.vertices[i3]->y, polySec.vertices[j4]->x, polySec.vertices[j4]->y)) {
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
                    cout << "(" << x << "; " << y << ")" << endl;
                    /*long double x = ((vertices[i1]->y - ((vertices[j2]->y - vertices[i1]->y) / (vertices[j2]->x - vertices[i1]->x))*vertices[i1]->x) -
                        (polySec.vertices[i3]->y - ((polySec.vertices[j4]->y - polySec.vertices[i3]->y) / (polySec.vertices[j4]->x - polySec.vertices[i3]->x))*polySec.vertices[i3]->x))
                        / (((polySec.vertices[j4]->y - polySec.vertices[i3]->y) / (polySec.vertices[j4]->x - polySec.vertices[i3]->x))
                            -((vertices[j2]->y - vertices[i1]->y) / (vertices[j2]->x - vertices[i1]->x)));
                    long double y = ((vertices[i1]->y - ((vertices[j2]->y - vertices[i1]->y) / (vertices[j2]->x - vertices[i1]->x))*vertices[i1]->x) -
                        (polySec.vertices[i3]->y - ((polySec.vertices[j4]->y - polySec.vertices[i3]->y) / (polySec.vertices[j4]->x - polySec.vertices[i3]->x))*polySec.vertices[i3]->x))
                        / (((polySec.vertices[j4]->y - polySec.vertices[i3]->y) / (polySec.vertices[j4]->x - polySec.vertices[i3]->x))
                            - ((vertices[j2]->y - vertices[i1]->y) / (vertices[j2]->x - vertices[i1]->x)))
                        * ((vertices[j2]->y - vertices[i1]->y) / (vertices[j2]->x - vertices[i1]->x)) + (vertices[i1]->y - ((vertices[j2]->y - vertices[i1]->y) / (vertices[j2]->x - vertices[i1]->x))*vertices[i1]->x);*/

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
        vector<shared_ptr<Point>> intersectVertices = pointIntersect(polySec);
        
        intersectionPolygon.vertices.insert(intersectionPolygon.vertices.end(), intersectVertices.begin(), intersectVertices.end());
        orderPointsClockwise(intersectionPolygon.vertices);

        return intersectionPolygon;
    }

    void orderPointsClockwise(vector<shared_ptr<Point>>& points) {
        //auto start = points[0];
        //for (const auto& vertex : points) {
        //    if (vertex->y < start->y || (vertex->y == start->y && vertex->x < start->x)) {
        //        start = vertex;
        //    }
        //}

        //// Сортируем оставшиеся точки по полярному углу относительно начальной точки
        //sort(points.begin(), points.end(), [start](const auto& a, const auto& b) {
        //    return atan2(a->y - start->y, a->x - start->x) < atan2(b->y - start->y, b->x - start->x);
        //    });
        auto start = points[0];
        for (const auto& vertex : points) {
            if (vertex->y < start->y || (vertex->y == start->y && vertex->x < start->x)) {
                start = vertex;
            }
        }

        // Сортируем оставшиеся точки по углу относительно начальной точки
        sort(points.begin(), points.end(), [start](const auto& a, const auto& b) {
            double angleA = atan2(a->y - start->y, a->x - start->x);
            double angleB = atan2(b->y - start->y, b->x - start->x);
            return angleA < angleB || (angleA == angleB && a->x < b->x);
            });
    }


    


    polygon unionPolygons(polygon polySec) {
        polygon unionPolygon(true);
        int count_1 = vertices.size();
        int count_2 = polySec.vertices.size();

        for (int i = 0; i < count_1; i++) {
            if (polySec.pointInside(vertices[i]->x, vertices[i]->y)==false) {
                unionPolygon.vertices.push_back(vertices[i]);
            }
        }
        for (int j = 0; j < count_2; j++) {
            if (pointInside(polySec.vertices[j]->x, polySec.vertices[j]->y)==false) {
                unionPolygon.vertices.push_back(polySec.vertices[j]);
            }
        }
        vector<shared_ptr<Point>> intersectVertices = pointIntersect(polySec);

        unionPolygon.vertices.insert(unionPolygon.vertices.end(), intersectVertices.begin(), intersectVertices.end());
        orderPointsClockwise(unionPolygon.vertices);

        return unionPolygon;
    }


  
    //void printPolygon() const {
    //    if (vertices.empty()) {
    //        cout << "Polygon is empty." << endl;
    //        return;
    //    }

    //    // Найти минимальные и максимальные координаты для определения размера поля вывода
    //    double minX = vertices[0]->x, minY = vertices[0]->y, maxX = vertices[0]->x, maxY = vertices[0]->y;
    //    for (const auto& vertex : vertices) {
    //        minX = (vertex->x < minX) ? vertex->x : minX;
    //        minY = (vertex->y < minY) ? vertex->y : minY;
    //        maxX = (vertex->x > maxX) ? vertex->x : maxX;
    //        maxY = (vertex->y > maxY) ? vertex->y : maxY;
    //    }

    //    // Создать двумерный массив для поля вывода
    //    vector<vector<char>> display(maxY - minY + 1, std::vector<char>(maxX - minX + 1, ' '));

    //    // Заполнить массив символами '*' на местах вершин многоугольника
    //    for (const auto& vertex : vertices) {
    //        display[static_cast<int>(vertex->y - minY)][static_cast<int>(vertex->x - minX)] = '*';
    //    }

    //    // Вывести массив на экран
    //    for (const auto& row : display) {
    //        for (char cell : row) {
    //            cout << cell;
    //        }
    //        cout << endl;
    //    }
    //}


    





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
    poly2.addVertex(1.5, 3);
    poly2.addVertex(3, 3);
    poly2.addVertex(3, 1);

    poly2.printPolygonVertex();


    poly.pointIntersect(poly2);*/



    /*polygon poly(true);
    poly.addVertex(0, 0);
    poly.addVertex(0, 4);
    poly.addVertex(4, 4);
    poly.addVertex(4, 0);

    poly.printPolygonVertex();

    polygon poly2(true);
    poly2.addVertex(1, 2);
    poly2.addVertex(7, 0);
    poly2.addVertex(7, 6);
    poly2.addVertex(3, 6);
    poly2.addVertex(3, 3);

    poly2.printPolygonVertex();*/

    polygon poly(true);
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

    poly2.printPolygonVertex();

    polygon polyInter(true);

    polyInter = poly.intersectionPolygons(poly2);
    polyInter.printPolygonVertex();

    polygon polyUni(true);

    polyUni = poly.unionPolygons(poly2);
    polyUni.printPolygonVertex();
    //poly.pointIntersect(poly2);


    /*cout << "Area: " << poly.calculateArea() << std::endl;
    cout << "Perimeter: " << poly.calculatePerimeter() << std::endl;*/



    /*cout << poly.pointInside(1, 1) << endl;
    cout << poly.pointInside(2, 2) << endl;
    cout << poly.pointInside(2, 0) << endl;
    cout << poly.pointInside(-1, 1) << endl;
    cout << poly.pointInside(3, 1) << endl;*/


    return 0;
}
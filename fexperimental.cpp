#include "fexperimental.h"
#include "qdir.h"
#include "qregularexpression.h"

FExperimental::FExperimental(QObject *parent)
    : QObject{parent}
{}

void FExperimental::ReadXY(const FCalculation &calculate)
{
    QFile file(path_exp_file);      //Открытие файла с экспериментальным спектром
    //Проверка открытия файла
    if (!file.open(QIODevice::ReadOnly)) {
        qWarning("Cannot open file for reading");
        exit(EXIT_FAILURE);
    }
    //Очистка старых значений координат
    expAbscissaX.clear();
    expOrdinateY.clear();
    QTextStream stream(&file);      //Поток данных из файла
    QStringList list;               //Список слов получаемых из строк
    //Заполнение новыми значениями
    while(!stream.atEnd()) {
        list = stream.readLine().split(QRegularExpression("\\s+"));
        expAbscissaX.push_back(list[0].toDouble());
        expOrdinateY.push_back(list[1].toDouble());
    }
    file.close();                   //Закрытие файла
    //Вычитание минимального значения из всего спектра
    double min_value_exp_y = *std::min_element(expOrdinateY.begin(), expOrdinateY.end());
    for (int i = 0; i < expOrdinateY.size(); ++i) {
        if (calculate.getSource() == FCalculation::HE_I) {
            expOrdinateY[i] -= (0.0 * min_value_exp_y);
        } else {
            expOrdinateY[i] -= (1.0 * min_value_exp_y);
        }
    }
    max_value_exp_y = *std::max_element(expOrdinateY.begin(), expOrdinateY.end());
}

const double &FExperimental::getMax_value_exp_y() const
{
    return max_value_exp_y;
}

void FExperimental::setMax_value_exp_y(double newMax_value_exp_y)
{
    max_value_exp_y = newMax_value_exp_y;
}

const QString &FExperimental::getPath_exp_file() const
{
    return path_exp_file;
}

void FExperimental::setPath_exp_file(const QString &newPath_exp_file)
{
    path_exp_file = newPath_exp_file;
}

const QList<double> &FExperimental::getExpOrdinateY() const
{
    return expOrdinateY;
}

void FExperimental::setExpOrdinateY(const QList<double> &newExpOrdinateY)
{
    expOrdinateY = newExpOrdinateY;
}

const QList<double> &FExperimental::getExpAbscissaX() const
{
    return expAbscissaX;
}

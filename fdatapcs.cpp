#include "fdatapcs.h"
#include "qdir.h"

FDataPCS::FDataPCS() {}

const QString FDataPCS::pcs_path = "./PhotoionizationCrossSection&Asymmetry.dat";

void FDataPCS::MiningPCS(const FOutput &output, QList<FElement *> &element)
{
    photoionizationCrossSection.clear();            //Очистка и перезаполнение вектора photoionizationCrossSection
    asymmetry.clear();
    photoionizationCrossSection.resize(4);
    asymmetry.resize(4);
    for (auto& sourse : photoionizationCrossSection) {
        sourse.resize(output.getUnique_atoms());
        for (auto& a : sourse) {
            a.resize(output.getKSpdfghi());
            for (auto& s : a)
                s.resize(19);
        }
    }
    for (auto& sourse : asymmetry) {
        sourse.resize(output.getUnique_atoms());
        for (auto& a : sourse) {
            a.resize(output.getKSpdfghi());
            for (auto& s : a)
                s.resize(19);
        }
    }

    QFile file(pcs_path);                           //Открытие файла с сечениями
    if (!file.open(QIODevice::ReadOnly)) {
        qWarning("Cannot open file for reading");   //Если файл не найден, то выводим предупреждение и завершаем выполнение программы
        exit(EXIT_FAILURE);
    }

    QTextStream stream(&file);
    QStringList list;   //Список слов получаемых из строк
    unsigned short a = 0, s = 0, n = 0;
    do {
        list = stream.readLine().split(u'\t');
        if (QString::compare(list[0], output.getUElements()[a], Qt::CaseInsensitive) == 0 || QString::compare(list[2], output.getUElements()[a], Qt::CaseInsensitive) == 0) {
            element[a]->setElementNumber(list[1].toUShort());
            element[a]->setElementName(list[2]);
            element[a]->setElementColor(QColor(list[3].toInt(), list[4].toInt(), list[5].toInt()));
            while (true) {
                list = stream.readLine().split(u'\t');
                if (list[0] == "asymmetry") {
                    break;
                } else {
                    n = list[0][0].digitValue() - 1;
                    switch (static_cast<Shells>(list[0][1].unicode())) {
                    case S:
                        s = 0;
                        break;
                    case P:
                        n -= 1;
                        s = 1;
                        break;
                    case D:
                        n -= 2;
                        s = 2;
                        break;
                    case F:
                        n -= 3;
                        s = 3;
                        break;
                    default:
                        break;
                    }
                    list.erase(list.cbegin());
                    list.erase(list.cend() - 1);
                    for (const auto &source : list)
                        photoionizationCrossSection[&source - &list[0]][a][s][n] +=  source.toDouble();
                }
            }
            while (true) {
                list = stream.readLine().split(u'\t');
                if (list[0] == '-') {
                    break;
                } else {
                    n = list[0][0].digitValue() - 1;
                    switch (static_cast<Shells>(list[0][1].unicode())) {
                    case S:
                        s = 0;
                        break;
                    case P:
                        n -= 1;
                        s = 1;
                        break;
                    case D:
                        n -= 2;
                        s = 2;
                        break;
                    case F:
                        n -= 3;
                        s = 3;
                        break;
                    default:
                        break;
                    }
                    list.erase(list.cbegin());
                    list.erase(list.cend() - 1);
                    for (const auto &source : list)
                        asymmetry[&source - &list[0]][a][s][n] +=  source.toDouble();
                }
            }
            a++;
        }
    } while (a < output.getUnique_atoms());

    QList<unsigned short> temp_element_number;
    for (int i = 0; i < output.getUnique_atoms(); i++)
        temp_element_number.push_back(element[i]->getElementNumber());
    for (int i = 0; i < output.getUnique_atoms(); i++) {    //Сортировка столбцов пузырьком
        for (int j = 0; j < output.getUnique_atoms() - 1; j++) {
            if (temp_element_number[j] < temp_element_number[j + 1]) {
                temp_element_number.swapItemsAt(j, j + 1);
                for (unsigned short source = 0; source < 4; source++) {
                    photoionizationCrossSection[source].swapItemsAt(j, j + 1);
                    asymmetry[source].swapItemsAt(j, j + 1);
                }
            }
        }
    }

    file.close();
}

const QList<QList<QList<QList<double> > > > &FDataPCS::getPhotoionizationCrossSection() const
{
    return photoionizationCrossSection;
}

const QList<QList<QList<QList<double> > > > &FDataPCS::getAsymmetry() const
{
    return asymmetry;
}

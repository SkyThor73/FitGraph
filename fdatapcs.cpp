#include "fdatapcs.h"
#include <QDir>
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QHash>
#include <QColor>

// -------------------------------------------------------------------------
// Расширенный статический кэш базы данных
// -------------------------------------------------------------------------
namespace {
struct PCSData {
    QList<double> sigma;
    QList<double> beta;
};

// Структура для хранения метаданных атома
struct ElementCache {
    unsigned short Z;
    QString name;
    QString colorHex;
    QHash<QString, PCSData> subshells;
};

QHash<QString, ElementCache> databaseCache;
bool isCacheLoaded = false;
}

const QString FDataPCS::pcs_path = "atomic_data.json";

FDataPCS::FDataPCS() {}

//Статический метод инициализации (Загружает JSON только один раз)
void FDataPCS::initializeDatabase()
{
    if (!isCacheLoaded) {
        QFile file(pcs_path);
        if (file.open(QIODevice::ReadOnly)) {
            QJsonDocument jsonDoc = QJsonDocument::fromJson(file.readAll());
            QJsonObject elementsObj = jsonDoc.object()["elements"].toObject();
            for (auto it = elementsObj.begin(); it != elementsObj.end(); ++it) {
                QString elemKey = it.key();
                QJsonObject elemObj = it.value().toObject();
                ElementCache cacheItem;
                cacheItem.Z = elemObj["Z"].toInt();
                cacheItem.name = elemKey;
                cacheItem.colorHex = elemObj["color"].toString();
                QJsonObject subshellsObj = elemObj["subshells"].toObject();
                for (auto shellIt = subshellsObj.begin(); shellIt != subshellsObj.end(); ++shellIt) {
                    QJsonObject shellObj = shellIt.value().toObject();
                    QJsonArray sigmaArr = shellObj["sigma"].toArray();
                    //Определение симметрии подоболочки
                    bool is_s_shell = shellIt.key().endsWith(u's');
                    QJsonArray betaArr;
                    //Извлечение массива beta только для p, d, f подоболочек
                    if (!is_s_shell)
                        betaArr = shellObj["beta"].toArray();
                    PCSData data;
                    const int target_indices[4] = {2, 4, 13, 14};
                    for (int idx : target_indices) {
                        data.sigma.push_back(idx < sigmaArr.size() ? sigmaArr[idx].toDouble() : 0.0);
                        //Аппаратное присвоение константы для s-электронов
                        if (is_s_shell) {
                            data.beta.push_back(2.0);
                        } else {
                            data.beta.push_back(idx < betaArr.size() ? betaArr[idx].toDouble() : 0.0);
                        }
                    }
                    cacheItem.subshells.insert(shellIt.key(), data);
                }
                databaseCache.insert(elemKey, cacheItem);
            }
            file.close();
            isCacheLoaded = true;
        }
    }
}

unsigned short FDataPCS::getAtomicNumber(const QString &symbol)
{
    initializeDatabase();
    auto it = databaseCache.find(symbol);
    if (it != databaseCache.end()) return it.value().Z;
    return 0; //Для неизвестных элементов (будут помещены в конец списка)
}

void FDataPCS::MiningPCS(const FOutput &output, QList<FElement *> &element)
{
    initializeDatabase(); //Гарантия загрузки кэша

    const unsigned short unique_atoms = output.getUnique_atoms();
    const unsigned short subshells = output.getKSpdfghi();
    const unsigned short max_n = 19;
    const unsigned short sources = 4;

    photoionizationCrossSection.assign(sources, QList<QList<QList<double>>>(unique_atoms, QList<QList<double>>(subshells, QList<double>(max_n, 0.0))));
    asymmetry.assign(sources, QList<QList<QList<double>>>(unique_atoms, QList<QList<double>>(subshells, QList<double>(max_n, 0.0))));

    const char shell_chars[] = "spdf";

    for (unsigned short a = 0; a < unique_atoms; ++a) {
        QString elemName = output.getUElements()[a];

        auto elemIt = databaseCache.find(elemName);
        if (elemIt == databaseCache.end()) continue;

        const ElementCache &eCache = elemIt.value();

        if (a < element.size() && element[a] != nullptr) {
            element[a]->setElementNumber(eCache.Z);
            element[a]->setElementName(eCache.name);
            element[a]->setElementColor(QColor(eCache.colorHex));
        }

        for (unsigned short s = 0; s < subshells; ++s) {
            for (unsigned short n = 0; n < max_n; ++n) {
                QString shellKey = QString::number(n + s + 1) + shell_chars[s];
                auto shellIt = eCache.subshells.find(shellKey);
                if (shellIt != eCache.subshells.end()) {
                    const PCSData &data = shellIt.value();
                    for (unsigned short source = 0; source < sources; ++source) {
                        if (source < data.sigma.size()) photoionizationCrossSection[source][a][s][n] = data.sigma[source];
                        if (source < data.beta.size()) asymmetry[source][a][s][n] = data.beta[source];
                    }
                }
            }
        }
    }
}

const QList<QList<QList<QList<double> > > > &FDataPCS::getPhotoionizationCrossSection() const
{
    return photoionizationCrossSection;
}

const QList<QList<QList<QList<double> > > > &FDataPCS::getAsymmetry() const
{
    return asymmetry;
}
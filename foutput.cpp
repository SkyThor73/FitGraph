#include "foutput.h"
#include "Core"
#include "Eigenvalues"
#include "qdir.h"
#include "qregularexpression.h"

FOutput::FOutput(QObject *parent)
    : QObject{parent}
{}

class PopulationAnalysis {
public:
    //Расчет полной матрицы заселенностей по Малликену
    //Возвращает матрицу, где элемент (mu, j) - вклад функции mu в орбиталь j
    static Eigen::MatrixXd calculateMulliken(const Eigen::MatrixXd& C, const Eigen::MatrixXd& S) {
        Eigen::MatrixXd SC = S * C;
        //Поэлементное умножение (Hadamard product)
        return C.cwiseProduct(SC);
    }

    //Расчет матрицы заселенностей по Лёвдину
    static Eigen::MatrixXd calculateLowdin(const Eigen::MatrixXd& C, const Eigen::MatrixXd& S) {
        //Диагонализация симметричной положительно определенной матрицы S
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);

        //S = U * D * U^T  => S^{1/2} = U * D^{1/2} * U^T
        Eigen::VectorXd D_sqrt = solver.eigenvalues().cwiseSqrt();
        Eigen::MatrixXd S_half = solver.eigenvectors() * D_sqrt.asDiagonal() * solver.eigenvectors().transpose();

        //Ортогонализация коэффициентов
        Eigen::MatrixXd C_ortho = S_half * C;

        //Возведение в квадрат каждого элемента
        return C_ortho.cwiseProduct(C_ortho);
    }
};

void FOutput::ReadData(const QString &path_file)
{
    QFile file(path_file);

    if (!file.open(QIODevice::ReadOnly)) {
        qWarning("Cannot open file for reading");
        exit(EXIT_FAILURE);
    }

    programm = NONE, runtype = false, mulliken = false, ecp = false, lowdin = false, sym = false; //Сброс ключей

    QTextStream stream(&file);
    QString str;

    do {                                 //Поиск программы в которой был получен выходной файл
        stream >> str;
        if (str == "GAMESS") {
            programm = GAMESS;
        } else if (str == "Firefly") {
            programm = FIREFLY;
        } else if (str == "Gaussian") {
            programm = GAUSSIAN;
        } else if (str == "*****************") {
            programm = ORCA;
        } else if (str == "CP2K") {
            programm = CP2K;
        }
    } while (programm == NONE);


    switch (programm) {
    case GAMESS:
    case FIREFLY:
        ReadingGamessFly(stream, str);
        break;
    case GAUSSIAN:
        ReadingGaussian(stream, str);
        break;
    case ORCA:
        ReadingOrca(stream, str);
        break;
    case CP2K:
        ReadingCP2K(stream, str);
        break;
    default:
        break;
    }

    file.close();                       //Закрытие файла

    uElements.clear();
    uElements = elements;
    std::sort(uElements.begin(), uElements.end());
    const auto last = std::unique(uElements.begin(), uElements.end());
    uElements.erase(last, uElements.cend());
    unique_atoms = uElements.size();
}

void FOutput::ReadingOrca(QTextStream &stream, QString &str)
{
    QStringList list;   //Список слов получаемых из строк ----- Orbital basis set information -----
    QString line;
    while (true) {
        if (stream.readLine() == "----- Orbital basis set information -----") {
            list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
            basis = list.last();
            break;
        }
    }
    //Поиск ключей
    while (true) {
        if (stream.readLine() == "                                       INPUT FILE") {
            do {
                str = stream.readLine();
                if (!sym) {
                    sym = str.contains("sym", Qt::CaseInsensitive);
                }
                if (!mulliken && str.contains("largeprint", Qt::CaseInsensitive)) {
                    mulliken = true;
                    lowdin = true;
                }
            } while (!str.contains("*", Qt::CaseInsensitive));
            break;
        }
    }
    elements.clear();
    while (true) {
        if (stream.readLine() == "CARTESIAN COORDINATES (ANGSTROEM)") {
            stream.readLine();
            list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
            while (!list.isEmpty()) {
                elements.push_back(list[0]);
                list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
            }
            atom_count = elements.size();
            break;
        }
    }
    //Параметры молекулы
    while (true) {
        str = stream.readLine();
        if (str == "ECP PARAMETER INFORMATION") {
            ecp = true;
            ecpAtom.clear();
            ecpAtom.resize(atom_count);
            ecpElectrons.clear();
            ecpElectrons.resize(atom_count);
            unsigned short temp_replacing_electrons = 0;
            while (true) {
                list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                if (list.size() > 1) {
                    if (list[0] == "ECP")
                        break;
                    if (list[0] == "Group")
                        temp_replacing_electrons = list[7].toUShort();
                    if (list[0] == "Atom") {
                        QString temp_atom_number = "";
                        for (auto &i : list[1]) {
                            if (i.isDigit())
                                temp_atom_number.push_back(i);
                        }
                        ecpElectrons[temp_atom_number.toUShort()] = temp_replacing_electrons;
                        ecpAtom[temp_atom_number.toUShort()] = true;
                    }
                }
            }
        }
        if (str == "General Settings:") {
            basis_functions = 0;
            do {
                list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                if (!list.isEmpty()) {
                    if (list.last() == "UHF" || list.last() == "ROHF") {
                        uhf = true;
                    } else if(list[2] == "Electrons") {
                        orbital_count = list.last().toUShort() / 2;
                    } else if (list[1] == "Dimension") {
                        basis_functions = list.last().toUShort();
                    }
                }
            } while (basis_functions == 0);
            break;
        }
    }
    energy.clear();
    while (true) {
        if (stream.readLine() == "ORBITAL ENERGIES") {
            for (unsigned short i = 0; i < 3; ++i)
                stream.readLine();
            do {
                list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                if (sym) {
                    energy.push_back(fabs(list[list.size() - 2].toDouble()));
                } else {
                    energy.push_back(fabs(list.last().toDouble()));
                }
            } while (energy.size() != orbital_count);
            break;
        }
    }
    eigenVectors.clear();                           //Очистка предыдущих значений массива эйгенвекторов
    mullikenVectors.clear();
    lowdinVector.clear();
    eigenVectors.resize(orbital_count);
    for (auto &o : eigenVectors) {
        o.resize(atom_count);
        for (auto &a : o)
            a.resize(kSpdfghi);
    }
    mullikenVectors = eigenVectors;
    lowdinVector = eigenVectors;
    symmetry.clear();
    symmetry.resize(orbital_count);
    while (true) {
        if (stream.readLine() == "MOLECULAR ORBITALS (RHF, ROHF)") {
            stream.readLine();
            unsigned short o = 0, a = 0, s = 0;                         //Начало заполнения массива эйгенвекторов
            do {
                str = "";
                for (unsigned short j = 0; j < basis_functions + 4; ++j) {
                    line = stream.readLine().trimmed();
                    QRegularExpression re(R"(^(\S+)\s+(\S+)\s+(-?\d+\.\d{6})(?:\s*)(-?\d+\.\d{6})(?:\s*)(-?\d+\.\d{6})(?:\s*)(-?\d+\.\d{6})(?:\s*)(-?\d+\.\d{6})(?:\s*)(-?\d+\.\d{6})\s*$)");
                    QRegularExpressionMatch m = re.match(line);
                    list.clear();
                    if (m.hasMatch()) {
                        for (int i = 1; i <= 8; ++i)
                            list << m.captured(i);
                        //qDebug() << list.join(" ");  // OK!
                    } else {
                        //qDebug() << "Parsing error:" << line;
                    }
                    if (j > 3) {
                        //Atoms
                        if (str.isEmpty()) {
                            str = list[0];
                        } else if (list[0] != str) {
                            ++a;
                            str = list[0];
                        }
                        //Shells
                        if (s != 0 && list[1][list[1].size() - 1] == 's') {
                            s = 0;
                        } else if (s != 1 && list[1][1] == 'p') {
                            s = 1;
                        } else if (s != 2 && list[1][1] == 'd') {
                            s = 2;
                        } else if (s != 3 && list[1][1] == 'f') {
                            s = 3;
                        } else if (s != 4 && list[1][1] == 'g') {
                            s = 4;
                        } else if (s != 5 && list[1][1] == 'h') {
                            s = 5;
                        }
                        //Filling
                        list.erase(list.cbegin(), list.cend() - 6); //Удаление первых двух слов
                        for (const auto &i : list) {
                            if (o < orbital_count) {
                                if (ecp && ecpAtom[a] && eigenVectors[o][a][s].isEmpty()) {
                                    switch (ecpElectrons[a]) {
                                    case 28:
                                        switch (s) {
                                        case 0:
                                            for (unsigned short i = 0; i < 3; ++i) {
                                                eigenVectors[o][a][s].push_back(0.0);
                                                mullikenVectors[o][a][s].push_back(0.0);
                                                lowdinVector[o][a][s].push_back(0.0);
                                            }
                                            break;
                                        case 1:
                                            for (unsigned short i = 0; i < 6; ++i) {
                                                eigenVectors[o][a][s].push_back(0.0);
                                                mullikenVectors[o][a][s].push_back(0.0);
                                                lowdinVector[o][a][s].push_back(0.0);
                                            }
                                            break;
                                        case 2:
                                            for (unsigned short i = 0; i < 5; ++i) {
                                                eigenVectors[o][a][s].push_back(0.0);
                                                mullikenVectors[o][a][s].push_back(0.0);
                                                lowdinVector[o][a][s].push_back(0.0);
                                            }
                                            break;
                                        case 3:             //Нет f
                                            break;
                                        default:
                                            break;
                                        }
                                        break;
                                    default:
                                        break;
                                    }
                                }
                                eigenVectors[o][a][s].push_back(pow(i.toDouble(), 2));  //Заполнение массива собственных значений
                                mullikenVectors[o][a][s].push_back(0.0);
                                lowdinVector[o][a][s].push_back(0.0);
                            }
                            ++o;
                        }
                        o -= 6;
                    }
                }
                o += 6;
                a = 0;
            } while (o < orbital_count);
            break;
        }
    }
    while (true) {
        if (stream.readLine() == "MULLIKEN ORBITAL POPULATIONS PER MO") {
            stream.readLine();
            stream.readLine();
            unsigned short o = 0, a = 0, s = 0, n = 0, i = 0, l = 0;                         //Начало заполнения массива эйгенвекторов
            QString atom_number = "", shell_number = "";
            do {
                for (unsigned short i = 0; i < 4; ++i)
                    stream.readLine();
                QList <QString> tempA;          //Вектор атомов для поиска повторяющихся элементов
                QList <QString> tempS;          //Вектор оболочек для поиска повторяющихся элементов
                while (true) {
                    list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                    if (list.isEmpty())
                        break;
                    tempA.push_back(list[0]);
                    //Atoms
                    atom_number = "";
                    for (auto &i : list[0]) {
                        if (i.isDigit())
                            atom_number.push_back(i);
                    }
                    a = atom_number.toUShort();
                    //Shells
                    shell_number = "";
                    if (list[1][list[1].size() - 1] == 's') {
                        tempS.push_back(list[1]);
                        for (auto &i : list[1]) {
                            if (i.isDigit()) {
                                shell_number.push_back(i);
                            } else {
                                break;
                            }
                        }
                    } else {
                        for (auto &i : list[1]) {
                            if (i.isDigit()) {
                                shell_number.push_back(i);
                                tempS.push_back(list[1][0]);
                            } else {
                                tempS[tempS.size() - 1].push_back(list[1][1]);
                                break;
                            }
                        }
                    }
                    n = shell_number.toUShort() - 1;
                    if (tempA.size() > 1 && tempS.size() > 1 && tempA[tempA.size() - 2] == tempA.last() && tempS[tempS.size() - 2] == tempS.last()) {
                        ++i;
                    } else {
                        i = 0;
                    }
                    if (list[1][list[1].size() - 1] == 's') {
                        s = 0;
                        l = n;
                    } else if (list[1][1] == 'p') {
                        s = 1;
                        l = 3 * n + i;
                    } else if (list[1][1] == 'd') {
                        s = 2;
                        l = 5 * n + i;
                    } else if (list[1][1] == 'f') {
                        s = 3;
                        l = 7 * n + i;
                    } else if (list[1][1] == 'g') {
                        s = 4;
                        l = 9 * n + i;
                    } else if (list[1][1] == 'h') {
                        s = 5;
                    }
                    //Filling
                    list.erase(list.cbegin(), list.cend() - 6); //Удаление первых двух слов
                    for (const auto &i : list) {
                        if (o < orbital_count)
                            mullikenVectors[o][a][s][l] += fabs(i.toDouble() * 0.01);  //Заполнение массива собственных значений
                        ++o;
                    }
                    o -= 6;
                }
                o += 6;
            } while (o < orbital_count);
            break;
        }
    }
    while (true) {
        if (stream.readLine() == "LOEWDIN ORBITAL POPULATIONS PER MO") {
            stream.readLine();
            stream.readLine();
            unsigned short o = 0, a = 0, s = 0, n = 0, i = 0, l = 0;                         //Начало заполнения массива эйгенвекторов
            QString atom_number = "", shell_number = "";
            do {
                for (unsigned short i = 0; i < 4; ++i)
                    stream.readLine();
                QList <QString> tempA;          //Вектор атомов для поиска повторяющихся элементов
                QList <QString> tempS;          //Вектор оболочек для поиска повторяющихся элементов
                while (true) {
                    list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                    if (list.isEmpty())
                        break;
                    tempA.push_back(list[0]);
                    //Atoms
                    atom_number = "";
                    for (auto &i : list[0]) {
                        if (i.isDigit())
                            atom_number.push_back(i);
                    }
                    a = atom_number.toUShort();
                    //Shells
                    shell_number = "";
                    if (list[1][list[1].size() - 1] == 's') {
                        tempS.push_back(list[1]);
                        for (auto &i : list[1]) {
                            if (i.isDigit()) {
                                shell_number.push_back(i);
                            } else {
                                break;
                            }
                        }
                    } else {
                        for (auto &i : list[1]) {
                            if (i.isDigit()) {
                                shell_number.push_back(i);
                                tempS.push_back(list[1][0]);
                            } else {
                                tempS[tempS.size() - 1].push_back(list[1][1]);
                                break;
                            }
                        }
                    }
                    n = shell_number.toUShort() - 1;
                    if (tempA.size() > 1 && tempS.size() > 1 && tempA[tempA.size() - 2] == tempA.last() && tempS[tempS.size() - 2] == tempS.last()) {
                        ++i;
                    } else {
                        i = 0;
                    }
                    if (list[1][list[1].size() - 1] == 's') {
                        s = 0;
                        l = n;
                    } else if (list[1][1] == 'p') {
                        s = 1;
                        l = 3 * n + i;
                    } else if (list[1][1] == 'd') {
                        s = 2;
                        l = 5 * n + i;
                    } else if (list[1][1] == 'f') {
                        s = 3;
                        l = 7 * n + i;
                    } else if (list[1][1] == 'g') {
                        s = 4;
                        l = 9 * n + i;
                    } else if (list[1][1] == 'h') {
                        s = 5;
                    }
                    //Filling
                    list.erase(list.cbegin(), list.cend() - 6); //Удаление первых двух слов
                    for (const auto &i : list) {
                        if (o < orbital_count)
                            lowdinVector[o][a][s][l] += i.toDouble() * 0.01;  //Заполнение массива собственных значений
                        ++o;
                    }
                    o -= 6;
                }
                o += 6;
            } while (o < orbital_count);
            break;
        }
    }
}

void FOutput::ReadingCP2K(QTextStream &stream, QString &str)
{
    static const QRegularExpression re(R"(^(?:MO\|\s+)?(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(-?\d+\.\d{6})(?:\s+(-?\d+\.\d{6}))?(?:\s+(-?\d+\.\d{6}))?(?:\s+(-?\d+\.\d{6}))?$)");
    bool found_basis = false, found_funcs = false, found_atoms = false, found_orbitals = false; //Флаги
    while (!stream.atEnd()) {
        str = stream.readLine();
        if (!found_basis && str.contains("GLOBAL| Basis set file name", Qt::CaseInsensitive)) { //Название базиса
            basis = str.section(u' ', -1, -1, QString::SectionSkipEmpty);
            found_basis = true;
        } else if (!found_funcs && str.contains("- Spherical basis functions:", Qt::CaseInsensitive)) { //Количество базисных функций
            basis_functions = str.section(u' ', -1, -1, QString::SectionSkipEmpty).toUShort();
            found_funcs = true;
        } else if (!found_atoms && str.contains("MODULE QUICKSTEP: ATOMIC COORDINATES IN ANGSTROM", Qt::CaseInsensitive)) { //Список химических элементов и количество атомов
            elements.clear();
            stream.readLine(); //Пропуск пустой строки
            stream.readLine(); //Пропуск строки заголовков
            QStringList list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
            while (!list.isEmpty()) {
                if (list.size() > 2)
                    elements.push_back(list[2]);
                list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
            }
            atom_count = elements.size();
            found_atoms = true;
        } else if (!found_orbitals && str.contains("Number of occupied orbitals:", Qt::CaseInsensitive)) { //Количество орбиталей
            orbital_count = str.section(u' ', -1, -1, QString::SectionSkipEmpty).toUShort();
            found_orbitals = true;
        }
        if (found_basis && found_funcs && found_atoms && found_orbitals) //Проверка полноты данных
            break;
    }
    if (!found_basis || !found_funcs || !found_atoms || !found_orbitals)
        throw std::runtime_error("CP2K parsing error: missing fundamental system properties in output."); //Завершение программы при сбое парсинга в случае отсутствия каких-либо данных
    //Универсальная лямбда-функция для трехмерного заполнения векторов
    auto allocate_3d = [](auto& container, unsigned short dim1, unsigned short dim2, unsigned short dim3) {
        container.clear();
        container.resize(dim1);
        for (auto& d1 : container) {
            d1.resize(dim2);
            for (auto& d2 : d1)
                d2.resize(dim3);
        }
    };
    allocate_3d(overlap, orbital_count, atom_count, kSpdfghi);
    allocate_3d(eigenVectors, orbital_count, atom_count, kSpdfghi);
    energy.clear();
    symmetry.clear();
    symmetry.resize(orbital_count);
    //Структура для быстрого маппинга индексов строк на атом и подуровень
    struct BasisMap {unsigned short atom_idx = 0, shell_idx = 0;};
    //Маска для заполнения заселённостей
    QList<BasisMap> basis_mask(basis_functions);
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(basis_functions, basis_functions); //Матрица перектрытия
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(basis_functions, orbital_count); //Матрица собственных векторов
    stream.seek(0); //Сброс указателя потока в начало файла
    while (!stream.atEnd()) {
        QString line = stream.readLine();
        if (line == " OVERLAP MATRIX") {
            unsigned short o = 0;
            while (o < basis_functions) {
                stream.readLine(); //Пропуск пустой строки
                stream.readLine(); //Пропуск строки заголовков
                const unsigned short block_start = o;
                for (unsigned short j = 0; j < basis_functions; ++j) {
                    line = stream.readLine().trimmed();
                    if (line == "") line = stream.readLine().trimmed(); //Пропуск пустых строк
                    QRegularExpressionMatch m = re.match(line);
                    if (m.hasMatch()) {
                        o = block_start;
                        for (int i = 5; i <= m.lastCapturedIndex(); ++i)
                            S(j, o++) = m.capturedView(i).toDouble(); //Заполнение матрицы перектрытия
                    }
                }
            }
        } else if (line == " MO| EIGENVALUES, OCCUPATION NUMBERS, AND SPHERICAL EIGENVECTORS") {
            unsigned short o = 0;
            while (o < orbital_count) {
                const unsigned short block_start = o;
                for (unsigned short j = 0; j < basis_functions + 3; ++j) {
                    line = stream.readLine().trimmed();
                    if (line == "MO|") line = stream.readLine().trimmed(); //Пропуск пустых строк
                    if (j == 1) { //Блок энергий
                        o = block_start; //Сброс счётчика орбиталей на первый столбец
                        for (auto token : QStringTokenizer{line, u' ', Qt::SkipEmptyParts} | std::views::drop(1)) {
                            if (o < orbital_count)
                                energy.push_back(token.toDouble() * -kHartree);
                            ++o;
                        }
                    } else if (j > 2) { //Блок коэффициентов
                        QRegularExpressionMatch m = re.match(line);
                        if (m.hasMatch()) {
                            unsigned short a = m.capturedView(2).toUShort() - 1; //Извлечение счетчика атома
                            unsigned short s = 0;
                            auto shell_type = m.capturedView(4);
                            const char c = shell_type.isEmpty() ? '\0' : shell_type[1].toLatin1();
                            switch (c) { //Извлечение счётчика подуровня
                            case 's': s = 0; break;
                            case 'p': s = 1; break;
                            case 'd': s = 2; break;
                            case 'f': s = 3; break;
                            case 'g': s = 4; break;
                            case 'h': s = 5; break;
                            default: continue;
                            }
                            if (o < 5 && (j - 3) < basis_functions)
                                basis_mask[j - 3] = {a, s}; //Заполнение маски
                            o = block_start; //Сброс счётчика орбиталей на первый столбец
                            for (int i = 5; i <= m.lastCapturedIndex(); ++i) {
                                if (o < orbital_count) {
                                    const double val = m.capturedView(i).toDouble();
                                    eigenVectors[o][a][s].push_back(val * val); //Заполнение квадратов собственных векторов
                                    C(j - 3, o) = val; //Заполнение собственных векторов
                                }
                                ++o;
                            }
                        }
                    }
                }
            }
        }
    }
    //Расчёт заселённостей по Малликену и Лёвдину
    Eigen::MatrixXd P_mull = PopulationAnalysis::calculateMulliken(C, S);
    Eigen::MatrixXd P_low = PopulationAnalysis::calculateLowdin(C, S);
    allocate_3d(mullikenVectors, orbital_count, atom_count, kSpdfghi);
    allocate_3d(lowdinVector, orbital_count, atom_count, kSpdfghi);
    for (unsigned short o = 0; o < orbital_count; ++o) {
        for (unsigned short j = 0; j < basis_functions; ++j) {
            const auto& map = basis_mask[j];
            mullikenVectors[o][map.atom_idx][map.shell_idx].push_back(P_mull(j, o));
            lowdinVector[o][map.atom_idx][map.shell_idx].push_back(P_low(j, o));
        }
    }
    mulliken = true; //Для активации кнопок
    lowdin = true; //Для активации кнопок
}

void FOutput::ReadingGamessFly(QTextStream &stream, QString &str)
{
    QStringList list;   //Список слов получаемых из строк
    //Поиск ключей
    while (true) {
        if (stream.readLine() == "            ECHO OF THE FIRST FEW INPUT CARDS -") {
            do {
                str = stream.readLine();
                if (!runtype)
                    runtype = str.contains("runtyp=optimize", Qt::CaseInsensitive);
                if (!cctyp) {
                    cctyp = str.contains("cctyp", Qt::CaseInsensitive);
                    mulliken = false;
                }
                if (!uhf)
                    uhf = str.contains("scftyp=uhf", Qt::CaseInsensitive);
                if (!cctyp && !mulliken)
                    mulliken = (str.contains("nprint=8", Qt::CaseInsensitive) || str.contains("nprint=9", Qt::CaseInsensitive));
                if (!lowdin)
                    lowdin = str.contains("nprint=9", Qt::CaseInsensitive);
                if (!ecp)
                    ecp = (str.contains("ecp=read", Qt::CaseInsensitive) || str.contains("pp=read", Qt::CaseInsensitive));
            } while (!str.contains("$end", Qt::CaseInsensitive));
            break;
        }
    }
    while (true) {
        list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
        if (list.size() > 0 && list[1] == "REQUESTS") {
            basis = list[4];
            break;
        } else if (list.size() > 0 && list[1] == "OPTIONS") {
            stream.readLine();
            list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
            list[0].erase(list[0].cbegin(), list[0].cbegin() + 7);
            basis = list[0];
            break;
        }
    }
    //Параметры молекулы
    while (true) {
        str = stream.readLine();
        if (str.contains("SHELLS")) {
            while (true) {
                list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                if (list[list.size() - 3] == "FUNCTIONS") {
                    basis_functions = list.last().toUShort();
                } else if(list[list.size() - 3] == "(ALPHA)") {
                    orbital_count = list.last().toUShort();
                } else if(list[list.size() - 3] == ")") {
                    orbital_count_beta = list.last().toUShort();
                } else if(list[list.size() - 3] == "ATOMS") {
                    atom_count = list.last().toUShort();
                    break;
                }
            }
            break;
        }
    }

    if (ecp) {                                      //ECP
        ecpAtom.clear();
        ecpAtom.resize(atom_count);
        ecpElectrons.clear();
        ecpElectrons.resize(atom_count);
        while (true) {
            if (stream.readLine() == "          ECP POTENTIALS") {
                while (true) {
                    list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                    if (list.size() > 1) {
                        if (list[0] == "THE")
                            break;
                        if (list[0] == "PARAMETERS") {
                            if (list[9] == "SAME") {
                                ecpElectrons[list[6].toUShort() - 1] = ecpElectrons[list[12].toUShort() - 1];
                                ecpAtom[list[6].toUShort() - 1] = ecpAtom[list[12].toUShort() - 1];
                            } else {
                                ecpElectrons[list[6].toUShort() - 1] = list[9].toUShort();
                                ecpAtom[list[6].toUShort() - 1] = true;
                            }
                        }
                    }
                }
                break;
            }
        }
        while (true) {
            if (uhf) {                              //Не отлажено
                if (stream.readLine() == " SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW.  ALPHA SET(S).") {
                    stream >> str;
                    orbital_count = str.toUShort();
                } else if (stream.readLine() == " SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW.   BETA SET(S).") {
                    stream >> str;
                    orbital_count_beta = str.toUShort();
                    break;
                }
            } else {
                if (stream.readLine() == " SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW.   BOTH SET(S).") {
                    stream >> str;
                    orbital_count = str.toUShort();
                    break;
                }
            }
        }
    }

    eigenVectors.clear();                           //Очистка предыдущих значений массива эйгенвекторов
    eigenVectors.resize(orbital_count);
    for (auto &o : eigenVectors) {
        o.resize(atom_count);
        for (auto &a : o)
            a.resize(kSpdfghi);
    }
    energy.clear();
    symmetry.clear();
    elements.clear();                                                   //Очистка предыдущих значений химических элементов
    while (true) {
        if (FoundItSimple(stream)) {
            unsigned short o = 0, a = 0, s = 0;                         //Начало заполнения массива эйгенвекторов
            do {
                stream.readLine();
                for (unsigned short j = 0; j < basis_functions + 3; ++j) {
                    list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                    if (j == 1) {                                       //В третьей строке информация об энергиях
                        for (const auto &i : list) {
                            if (o < orbital_count)                      //Сбор данных с заполненых орбиталей
                                energy.push_back(i.toDouble() * -kHartree); //Заполнение вектора энергий
                            o++;
                        }
                        o -= 5;
                    } else if (j == 2) {                                       //В третьей строке информация об энергиях
                        for (const auto &i : list) {
                            if (o < orbital_count)                      //Сбор данных по симметрии МО
                                symmetry.push_back(i);
                            o++;
                        }
                        o -= 5;
                    } else if (j > 2) {
                        //Atoms
                        if (list[1].size() == 4) {  //Если название атома из двух символов и атомов больше 10
                            list.insert(2, list[1]);
                            list[1].erase(list[1].cbegin() + 2, list[1].cend());
                            list[2].erase(list[2].cbegin(), list[2].cend() - 2);
                        } else if (list[1].size() == 5) {
                            list.insert(2, list[1]);        //Копируем значение в соседнее положение
                            list[1].erase(list[1].cbegin() + 2, list[1].cend());
                            list[2].erase(list[2].cbegin(), list[2].cend() - 3);
                        }
                        if ((a + 1) == (list[2].toUShort() - 1) || elements.size() == 0) {
                            a = list[2].toUShort() - 1;
                            s = 0;
                            //Elements
                            if (o == 0)
                                elements.push_back(list[1]);
                        }
                        //Shells
                        if (s != 0 && list[3] == "S") {
                            s = 0;
                        } else if (s != 1 && list[3] == "X") {
                            s = 1;
                        } else if (s != 2 && list[3] == "XX") {
                            s = 2;
                        } else if (s != 3 && list[3] == "XXX") {
                            s = 3;
                        } else if (s != 4 && list[2][list[2].size() - 1] == 'X') {
                            s = 4;
                        } else if (s != 5 && list[3] == "XXXXX") {
                            s = 5;
                        }
                        //Filling
                        list.erase(list.cbegin(), list.cend() - 5); //Удаление первых двух слов
                        for (const auto &i : list) {
                            if (o < orbital_count) {
                                if (ecp && ecpAtom[a] && eigenVectors[o][a][s].isEmpty()) {
                                    switch (ecpElectrons[a]) {
                                    case 28:
                                        switch (s) {
                                        case 0:
                                            for (unsigned short i = 0; i < 3; ++i)
                                                eigenVectors[o][a][s].push_back(0.0);
                                            break;
                                        case 1:
                                            for (unsigned short i = 0; i < 6; ++i)
                                                eigenVectors[o][a][s].push_back(0.0);
                                            break;
                                        case 2:
                                            for (unsigned short i = 0; i < 6; ++i)
                                                eigenVectors[o][a][s].push_back(0.0);
                                            break;
                                        case 3:             //Нет f
                                            break;
                                        default:
                                            break;
                                        }
                                        break;
                                    default:
                                        break;
                                    }
                                }
                                eigenVectors[o][a][s].push_back(pow(i.toDouble(), 2));  //Заполнение массива собственных значений
                            }
                            o++;
                        }
                        o -= 5;
                    }
                }
                o += 5;
                a = 0;
            } while (o < orbital_count);
            break;
        }
    }

    if (uhf) {                                                 //uhf
        eigenVectorsBeta.clear();
        eigenVectorsBeta.resize(orbital_count_beta);
        for (auto &o : eigenVectorsBeta) {
            o.resize(atom_count);
            for (auto &a : o)
                a.resize(kSpdfghi);
        }
        energyBeta.clear();
        symmetryBeta.clear();
        while (true) {
            if (FoundItSimple(stream)) {
                unsigned short o = 0, a = 0, s = 0;                         //Начало заполнения массива эйгенвекторов
                do {
                    stream.readLine();
                    for (unsigned short j = 0; j < basis_functions + 3; ++j) {
                        list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                        if (j == 1) {                                       //В третьей строке информация об энергиях
                            for (const auto &i : list) {
                                if (o < orbital_count_beta)                      //Сбор данных с заполненых орбиталей
                                    energyBeta.push_back(i.toDouble() * -kHartree); //Заполнение вектора энергий
                                o++;
                            }
                            o -= 5;
                        } else if (j == 2) {                                       //В третьей строке информация об энергиях
                            for (const auto &i : list) {
                                if (o < orbital_count_beta)                      //Сбор данных по симметрии МО
                                    symmetryBeta.push_back(i);
                                o++;
                            }
                            o -= 5;
                        } else if (j > 2) {
                            //Atoms
                            if (list[1].size() == 4) {  //Если название атома из двух символов и атомов больше 10
                                list.insert(2, list[1]);
                                list[1].erase(list[1].cbegin() + 2, list[1].cend());
                                list[2].erase(list[2].cbegin(), list[2].cend() - 2);
                            }
                            if ((a + 1) == (list[2].toUShort() - 1)) {
                                a = list[2].toUShort() - 1;
                                s = 0;
                            }
                            //Shells
                            if (s != 0 && list[3] == "S") {
                                s = 0;
                            } else if (s != 1 && list[3] == "X") {
                                s = 1;
                            } else if (s != 2 && list[3] == "XX") {
                                s = 2;
                            } else if (s != 3 && list[3] == "XXX") {
                                s = 3;
                            } else if (s != 4 && list[2][list[2].size() - 1] == 'X') {
                                s = 4;
                            } else if (s != 5 && list[3] == "XXXXX") {
                                s = 5;
                            }
                            //Filling
                            list.erase(list.cbegin(), list.cend() - 5); //Удаление первых двух слов
                            for (const auto &i : list) {
                                if (o < orbital_count_beta)
                                    eigenVectorsBeta[o][a][s].push_back(pow(i.toDouble(), 2));  //Заполнение массива собственных значений
                                o++;
                            }
                            o -= 5;
                        }
                    }
                    o += 5;
                    a = 0;
                } while (o < orbital_count_beta);
                break;
            }
        }
    }

    if (lowdin) {                                                //NPRINT=9 BEGIN
        lowdinVector.clear();                                //Очистка предыдущих значений массива малликеновских векторов
        lowdinVector.resize(orbital_count);
        for (auto &o : lowdinVector) {
            o.resize(atom_count);
            for (auto &a : o)
                a.resize(kSpdfghi);
        }
        while (true) {
            str = stream.readLine();
            if (FoundItLowdin(str, stream)) {
                unsigned short o = 0, a = 0, s = 0;                         //Начало заполнения массива
                do {
                    for(unsigned short i = 0; i < 5; ++i)
                        stream.readLine();
                    for (unsigned short j = 0; j < basis_functions; ++j) {
                        list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                        //Atoms
                        if (list[1].size() == 4) {  //Если название атома из двух символов и атомов больше 10
                            list.insert(2, list[1]);
                            list[1].erase(list[1].cbegin() + 2, list[1].cend());
                            list[2].erase(list[2].cbegin(), list[2].cend() - 2);
                        }
                        if ((a + 1) == (list[2].toUShort() - 1)) {
                            a = list[2].toUShort() - 1;
                            s = 0;
                        }
                        //Shells
                        if (s != 0 && list[3] == "S") {
                            s = 0;
                        } else if (s != 1 && list[3] == "X") {
                            s = 1;
                        } else if (s != 2 && list[3] == "XX") {
                            s = 2;
                        } else if (s != 3 && list[3] == "XXX") {
                            s = 3;
                        } else if (s != 4 && list[2][list[2].size() - 1] == 'X') {
                            s = 4;
                        } else if (s != 5 && list[3] == "XXXXX") {
                            s = 5;
                        }
                        //Filling
                        list.erase(list.cbegin(), list.cend() - 5); //Удаление первых двух слов
                        for (const auto &i : list) {
                            if (o < orbital_count) {
                                if (ecp && ecpAtom[a] && lowdinVector[o][a][s].isEmpty()) {
                                    switch (ecpElectrons[a]) {
                                    case 28:
                                        switch (s) {
                                        case 0:
                                            for (unsigned short i = 0; i < 3; ++i)
                                                lowdinVector[o][a][s].push_back(0.0);
                                            break;
                                        case 1:
                                            for (unsigned short i = 0; i < 6; ++i)
                                                lowdinVector[o][a][s].push_back(0.0);
                                            break;
                                        case 2:
                                            for (unsigned short i = 0; i < 6; ++i)
                                                lowdinVector[o][a][s].push_back(0.0);
                                            break;
                                        case 3:             //Не отлажено f
                                            break;
                                        default:
                                            break;
                                        }
                                        break;
                                    default:
                                        break;
                                    }
                                }
                                lowdinVector[o][a][s].push_back(fabs(i.toDouble()));  //Заполнение массива собственных значений
                            }
                            o++;
                        }
                        o -= 5;
                    }
                    o += 5;
                    a = 0;
                } while (o < orbital_count);
                break;
            }
        }
        if (uhf) {
            lowdinVectorBeta.clear();                                //Очистка предыдущих значений массива малликеновских векторов
            lowdinVectorBeta.resize(orbital_count_beta);
            for (auto &o : lowdinVectorBeta) {
                o.resize(atom_count);
                for (auto &a : o)
                    a.resize(kSpdfghi);
            }
            while (true) {
                str = stream.readLine();
                if (FoundItLowdin(str, stream)) {
                    unsigned short o = 0, a = 0, s = 0;                         //Начало заполнения массива
                    do {
                        for(unsigned short i = 0; i < 5; ++i)
                            stream.readLine();
                        for (unsigned short j = 0; j < basis_functions; ++j) {
                            list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                            //Atoms
                            if (list[1].size() == 4) {  //Если название атома из двух символов и атомов больше 10
                                list.insert(2, list[1]);
                                list[1].erase(list[1].cbegin() + 2, list[1].cend());
                                list[2].erase(list[2].cbegin(), list[2].cend() - 2);
                            }
                            if ((a + 1) == (list[2].toUShort() - 1)) {
                                a = list[2].toUShort() - 1;
                                s = 0;
                            }
                            //Shells
                            if (s != 0 && list[3] == "S") {
                                s = 0;
                            } else if (s != 1 && list[3] == "X") {
                                s = 1;
                            } else if (s != 2 && list[3] == "XX") {
                                s = 2;
                            } else if (s != 3 && list[3] == "XXX") {
                                s = 3;
                            } else if (s != 4 && list[2][list[2].size() - 1] == 'X') {
                                s = 4;
                            } else if (s != 5 && list[3] == "XXXXX") {
                                s = 5;
                            }
                            //Filling
                            list.erase(list.cbegin(), list.cend() - 5); //Удаление первых двух слов
                            for (const auto &i : list) {
                                if (o < orbital_count_beta)
                                    lowdinVectorBeta[o][a][s].push_back(fabs(i.toDouble()));  //Заполнение массива собственных значений
                                o++;
                            }
                            o -= 5;
                        }
                        o += 5;
                        a = 0;
                    } while (o < orbital_count_beta);
                    break;
                }
            }
        }
    }                               //NPRINT=9 END

    if (mulliken) {                                             //NPRINT=8 BEGIN
        QList<QList<double>> tempMullikenPopulations, tempMullikenPopulationsBeta;
        tempMullikenPopulations.resize(orbital_count);
        tempMullikenPopulationsBeta.resize(orbital_count_beta);

        while (true) {
            if (FoundItMulliken(stream)) { //Поиск строки
                unsigned short o = 0;
                do {
                    for(unsigned short i = 0; i < 5; ++i)
                        stream.readLine();
                    for (int b = 0; b < basis_functions; ++b) {
                        list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                        list.pop_front();
                        for (const auto &i : list) {
                            tempMullikenPopulations[o].push_back(fabs(i.toDouble()));
                            o++;
                        }
                        o -= list.size();
                    }
                    o += 5;
                } while (o < orbital_count);
                break;
            }
        }

        if (uhf) {                                                      //uhf
            while (true) {
                if (FoundItMulliken(stream)) { //Поиск строки
                    unsigned short o = 0;
                    do {
                        for(unsigned short i = 0; i < 5; ++i)
                            stream.readLine();
                        for (int b = 0; b < basis_functions; ++b) {
                            list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                            list.pop_front();
                            for (const auto &i : list) {
                                tempMullikenPopulationsBeta[o].push_back(fabs(i.toDouble()));
                                o++;
                            }
                            o -= list.size();
                        }
                        o += 5;
                    } while (o < orbital_count_beta);
                    break;
                }
            }
        }

        mullikenVectors.clear();                                //Очистка предыдущих значений массива малликеновских векторов
        mullikenVectors.resize(orbital_count);
        for (auto &o : mullikenVectors) {
            o.resize(atom_count);
            for (auto &a : o)
                a.resize(kSpdfghi);
        }
        if (uhf) {
            mullikenVectorsBeta.clear();                                //Очистка предыдущих значений массива малликеновских векторов бета
            mullikenVectorsBeta.resize(orbital_count_beta);
            for (auto &o : mullikenVectorsBeta) {
                o.resize(atom_count);
                for (auto &a : o)
                    a.resize(kSpdfghi);
            }
        }

        while (true) {
            if (stream.readLine() == "               ----- POPULATIONS IN EACH AO -----") {
                stream.readLine();
                unsigned short a = 0, s = 0;
                for (unsigned short b = 0; b < basis_functions; ++b) {
                    list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                    //Atoms
                    if (list[1].size() == 4) {  //Если название атома из двух символов и атомов больше 10
                        list.insert(2, list[1]);
                        list[1].erase(list[1].cbegin() + 2, list[1].cend());
                        list[2].erase(list[2].cbegin(), list[2].cend() - 2);
                    }
                    if ((a + 1) == (list[2].toUShort() - 1)) {
                        a = list[2].toUShort() - 1;
                        s = 0;
                    }
                    //Shells
                    if (s != 0 && list[3] == "S") {
                        s = 0;
                    } else if (s != 1 && list[3] == "X") {
                        s = 1;
                    } else if (s != 2 && list[3] == "XX") {
                        s = 2;
                    } else if (s != 3 && list[3] == "XXX") {
                        s = 3;
                    } else if (s != 4 && list[2][list[2].size() - 1] == 'X') {
                        s = 4;
                    } else if (s != 5 && list[3] == "XXXXX") {
                        s = 5;
                    }
                    //Filling
                    for (unsigned short o = 0; o < orbital_count; ++o) {
                        if (ecp && ecpAtom[a] && mullikenVectors[o][a][s].isEmpty()) {
                            switch (ecpElectrons[a]) {
                            case 28:
                                switch (s) {
                                case 0:
                                    for (unsigned short i = 0; i < 3; ++i)
                                        mullikenVectors[o][a][s].push_back(0.0);
                                    break;
                                case 1:
                                    for (unsigned short i = 0; i < 6; ++i)
                                        mullikenVectors[o][a][s].push_back(0.0);
                                    break;
                                case 2:
                                    for (unsigned short i = 0; i < 6; ++i)
                                        mullikenVectors[o][a][s].push_back(0.0);
                                    break;
                                case 3:             //Не отлажено f
                                    break;
                                default:
                                    break;
                                }
                                break;
                            default:
                                break;
                            }
                        }
                        mullikenVectors[o][a][s].push_back(tempMullikenPopulations[o][b] * 0.5);  //Заполнение малликеновских векторов
                    }
                    if (uhf) {
                        for (unsigned short o = 0; o < orbital_count_beta; ++o)
                            mullikenVectorsBeta[o][a][s].push_back(tempMullikenPopulationsBeta[o][b] * 0.5);  //Заполнение малликеновских векторов
                    }
                }
                break;
            }
        }
    }                               //NPRINT=8 END
}

void FOutput::ReadingGaussian(QTextStream &stream, QString &str)
{
    QStringList list;   //Список слов получаемых из строк
    while (true) {
        list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);    //Заполнение списка словами разделёнными пробелами
        if (list.size() > 2 && list[0] == "Standard") {    //Поиск нужной строки
            basis = list[2];
            break;
        }
    }
    while (true) {
        list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);    //Заполнение списка словами разделёнными пробелами
        if (list.size() > 3 && list[1] == "alpha") {    //Поиск нужной строки
            orbital_count = list[0].toUShort();         //Присвоение количества орбиталей альфа и бета
            orbital_count_beta = list[3].toUShort();
            break;
        }
    }
    while (true) {
        stream >> str;
        if (str == "NAtoms=") {     //Поиск количества атомов
            stream >> str;
            atom_count = str.toUShort();
            break;
        }
    }
    while (true) {
        stream >> str;
        if (str == "NBasis=") {     //Поиск количества базисных функций
            stream >> str;
            basis_functions = str.toUShort();
            break;
        }
    }

    eigenVectors.clear();                           //Очистка предыдущих значений массива эйгенвекторов
    eigenVectors.resize(orbital_count);
    for (auto &o : eigenVectors) {
        o.resize(atom_count);
        for (auto &a : o) {
            a.resize(kSpdfghi);
        }
    }
    simpleMu.clear();
    simpleMu.resize(orbital_count);
    for (auto &o : simpleMu) {
        o.resize(basis_functions);
    }
    energy.clear();
    elements.clear();                                                   //Очистка предыдущих значений химических элементов
    while (true) {
        str = stream.readLine();
        if (str == "     Molecular Orbital Coefficients:") {
            unsigned short o = 0, a = 0, s = 0;                         //Начало заполнения массива эйгенвекторов
            do {
                for (unsigned short var = 0; var < basis_functions + 3; ++var) {
                    list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                    if (var == 2) {                                     //В третьей строке информация об энергиях
                        list.erase(list.cbegin(), list.cbegin()+2);     //Удаление первых двух слов
                        for (const auto &i : list) {
                            if (o < orbital_count)                      //Сбор данных с заполненых орбиталей
                                energy.push_back(i.toDouble() * -kHartree); //Заполнение вектора энергий
                            o++;
                        }
                        o -= 5;
                    } else if (var > 2) {
                        if (list.size() == 7 || list.size() == 8) {     //Строки с Eigenvalues
                            if (list[1][list[1].size() - 1] == 'X') {   //Заменить это на switch добавить enum shell
                                s = 1;
                            } else if (list[1][list[1].size() - 1] == 'D') {
                                s = 2;
                                list.erase(list.cbegin());
                            } else if (list[1][list[1].size() - 1] == 'F') {
                                s = 3;
                                list.erase(list.cbegin());
                            }
                            list.erase(list.cbegin(), list.cbegin()+2); //Удаление первых двух слов
                        } else {                                        //Строки с атомами и Eigenvalues
                            a = list[1].toUShort() - 1;                 //Порядковый номер атома (только в этой строке)
                            if (o == 0)
                                elements.push_back(list[2]);            //Сокращённое название элемента
                            s = 0;                                      //Потому что эта строка всегда содержит 1S
                            list.erase(list.cbegin(), list.cbegin()+4); //Удаление первых четырёх слов
                        }
                        for (const auto &i : list) {
                            if (o < orbital_count) {
                                eigenVectors[o][a][s].push_back(pow(i.toDouble(), 2));  //Заполнение массива собственных значений
                                simpleMu[o][var - 3] = i.toDouble();
                            }
                            o++;
                        }
                        o -= 5;
                    }
                }
                o += 5;
            } while (o < orbital_count);
            break;
        }
    }
    mullikenMu.clear();
    mullikenMu.resize(orbital_count);
    for (auto &o : mullikenMu) {
        o.resize(basis_functions);
    }
    mulliken = true;
    overlapMatrix.clear();
    overlapMatrix.resize(basis_functions);
    while (true) {
        str = stream.readLine();
        if (str == "    Full Mulliken population analysis:") {
            do {
                stream.readLine();
                for (unsigned short var = 0; var < basis_functions; ++var) {
                    list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                    if (var == 0) {
                        var = list[0].toUShort() - 1;
                    }
                    if ((var - overlapMatrix[var].size()) / 5 == 0) {
                        list.erase(list.cbegin(), list.cend() - (var - overlapMatrix[var].size() + 1));
                        for (const auto &i : list) {
                            overlapMatrix[var].push_back(i.toDouble());
                        }
                    } else {
                        list.erase(list.cbegin(), list.cend() - 5);
                        for (const auto &i : list) {
                            overlapMatrix[var].push_back(i.toDouble());
                        }
                    }
                }
            } while (overlapMatrix[basis_functions - 1].size() < basis_functions);
            break;
        }
    }
    for (unsigned short o = 0; o < orbital_count; ++o) {
        for (unsigned short b = 0; b < basis_functions; ++b) {
            double temp_sum = 0;
            for (unsigned short nu = 0; nu < basis_functions; ++nu) {
                if (b < nu) {
                    temp_sum += (simpleMu[o][nu] * overlapMatrix[nu][b]);
                } else if (b > nu) {
                    temp_sum += (simpleMu[o][nu] * overlapMatrix[b][nu]);
                }
            }
            mullikenMu[o][b] = pow(simpleMu[o][b], 2) + simpleMu[o][b] * temp_sum;
        }
    }
    mullikenVectors.clear();                                //Очистка предыдущих значений массива малликеновских векторов
    mullikenVectors.resize(orbital_count);
    for (auto &o : mullikenVectors) {
        o.resize(atom_count);
        for (auto &a : o) {
            a.resize(kSpdfghi);
        }
    }
    while (true) {
        str = stream.readLine();
        if (str == "     Gross orbital populations:") {
            stream.readLine();
            unsigned short a = 0, s = 0;
            for (int b = 0; b < basis_functions; ++b) {
                list = stream.readLine().split(u' ',  Qt::SkipEmptyParts);
                if (list.size() == 5) {
                    a = list[1].toUShort() - 1;
                    s = 0;
                } else if (list[1][list[1].size() - 1] == 'X') {
                    s = 1;
                } else if (list[1][list[1].size() - 1] == 'D') {
                    s = 2;
                } else if (list[1][list[1].size() - 1] == 'F') {
                    s = 3;
                }
                for (unsigned short o = 0; o < orbital_count; o++) {
                    mullikenVectors[o][a][s].push_back(mullikenMu[o][b]);   //Заполнение малликеновских векторов
                }
            }
            break;
        }
    }
}

const FOutput::Programms &FOutput::getProgramm() const
{
    return programm;
}

const bool &FOutput::getUhf() const
{
    return uhf;
}

const bool &FOutput::getMulliken() const
{
    return mulliken;
}

const bool &FOutput::getLowdin() const
{
    return lowdin;
}

const unsigned short &FOutput::getKSpdfghi() const
{
    return kSpdfghi;
}

const unsigned short &FOutput::getOrbital_count() const
{
    return orbital_count;
}

const unsigned short &FOutput::getOrbital_count_beta() const
{
    return orbital_count_beta;
}

const unsigned short &FOutput::getAtom_count() const
{
    return atom_count;
}

const unsigned short &FOutput::getUnique_atoms() const
{
    return unique_atoms;
}

const QString &FOutput::getBasis() const
{
    return basis;
}

const QList<double> &FOutput::getEnergy() const
{
    return energy;
}

const QList<double> &FOutput::getEnergyBeta() const
{
    return energyBeta;
}

const QList<QString> &FOutput::getElements() const
{
    return elements;
}

const QList<QString> &FOutput::getUElements() const
{
    return uElements;
}

const QList<QString> &FOutput::getSymmetry() const
{
    return symmetry;
}

const QList<QString> &FOutput::getSymmetryBeta() const
{
    return symmetryBeta;
}

const QList<QList<QList<QList<double> > > > &FOutput::getEigenVectors() const
{
    return eigenVectors;
}

const QList<QList<QList<QList<double> > > > &FOutput::getEigenVectorsBeta() const
{
    return eigenVectorsBeta;
}

const QList<QList<QList<QList<double> > > > &FOutput::getMullikenVectors() const
{
    return mullikenVectors;
}

const QList<QList<QList<QList<double> > > > &FOutput::getMullikenVectorsBeta() const
{
    return mullikenVectorsBeta;
}

const QList<QList<QList<QList<double> > > > &FOutput::getLowdinVector() const
{
    return lowdinVector;
}

const QList<QList<QList<QList<double> > > > &FOutput::getLowdinVectorBeta() const
{
    return lowdinVectorBeta;
}

void FOutput::SwitchProg(const QString &str, unsigned short &i, unsigned short &a, unsigned short &s)
{
    switch (programm) {
    case GAMESS:
        if ((str[0] != 'X' || str[0] != 'Y' || str[0] != 'Z') && (str[str.size() - 1] == 'X' || str[str.size() - 1] == 'Y' || str[str.size() - 1] == 'Z')) {
            QString str_number;
            unsigned short x = 0;
            do {
                str_number.push_back(str[x]);
                x++;
            } while (str[x] != 'X' && str[x] != 'Y' && str[x] != 'Z');
            a = str_number.toUShort();
            s = str.size() - str_number.size();
            a--;
            i++;
        } else {
            a = str.toUShort();
            a--;
        }
        break;
    case FIREFLY:
        if (str[0] == 'X' || str[0] == 'Y' || str[0] == 'Z') {
            s = str.length();
            i++;
        } else {
            a = str.toUShort();
            a--;
        }
        break;
    default:
        break;
    }
}

bool FOutput::SwitchProg(const QString &str)
{
    switch (programm) {
    case GAMESS:
    case FIREFLY:
        if (str == "          MULLIKEN AND LOWDIN POPULATION ANALYSES") {
            return true;
        } else if (str == "         BETA ORBITALS") {
            return true;
        } else {
            return false;
        }
        break;
    default:
        return false;
        break;
    }
}

bool FOutput::FoundItSimple(QTextStream &stream)
{
    if (!runtype) {                             //runtype = ENERGY
        if (stream.readLine() == "          ------------") {
            stream.readLine();
            stream.readLine();
            return true;
        } else {
            return false;
        }
    } else {                                    //runtype = OPTIMIZE
        if (!uhf) {                                  //RHF
            if (stream.readLine() == "          MOLECULAR ORBITALS") {
                if (programm == FIREFLY) {
                    for (unsigned short i = 0; i < 6; i++)
                        stream.readLine();
                } else {
                    stream.readLine();
                }
                return true;
            } else {
                return false;
            }
        } else {                                //uhf
            if (stream.readLine() == "          **** ALPHA SET **** ") {
                if (programm == FIREFLY) {
                    for (unsigned short i = 0; i < 5; i++)
                        stream.readLine();
                }
                return true;
            } else if (stream.readLine() == "          **** BETA SET ****") {
                return true;
            } else {
                return false;
            }
        }
    }
}

bool FOutput::FoundItMulliken(QTextStream &stream)
{
    if (!runtype) {                             //runtype = ENERGY
        if (!uhf) {
            if (stream.readLine() == "          ---------------------------------------") {
                for(unsigned short i = 0; i < 4; ++i)
                    stream.readLine();
                return true;
            } else {
                return false;
            }
        } else {
            if (stream.readLine() == "     AO MULLIKEN POPULATIONS IN EACH MOLECULAR ORBITAL") {
                stream.readLine();
                return true;
            } else {
                return false;
            }
        }
    } else {                                    //runtype = OPTIMIZE
        if (!uhf) {                                  //RHF
            if (stream.readLine() == "          MOLECULAR ORBITALS") {
                if (programm == FIREFLY) {
                    for (unsigned short i = 0; i < 6; i++)
                        stream.readLine();
                } else {
                    stream.readLine();
                }
                return true;
            } else {
                return false;
            }
        } else {                                //uhf
            if (stream.readLine() == "          **** ALPHA SET **** ") {
                if (programm == FIREFLY) {
                    for (unsigned short i = 0; i < 5; i++)
                        stream.readLine();
                }
                return true;
            } else if (stream.readLine() == "          **** BETA SET ****") {
                return true;
            } else {
                return false;
            }
        }
    }
}

bool FOutput::FoundItLowdin(const QString &str, QTextStream &stream)
{
    if (!runtype) {                             //runtype = ENERGY
        if (!uhf) {
            if (str == "          -----------------------------------------" || str == "          ----------------------------------------------------------") {
                stream.readLine();
                stream.readLine();
                return true;
            } else {
                return false;
            }
        } else {
            if (str == "          ----------------------------------------------------------" || str == "          ---------------------------------------------------------") {
                stream.readLine();
                stream.readLine();
                return true;
            } else {
                return false;
            }
        }
    } else {                                    //runtype = OPTIMIZE
        if (!uhf) {                                  //RHF
            if (stream.readLine() == "          MOLECULAR ORBITALS") {
                if (programm == FIREFLY) {
                    for (unsigned short i = 0; i < 6; i++)
                        stream.readLine();
                } else {
                    stream.readLine();
                }
                return true;
            } else {
                return false;
            }
        } else {                                //uhf
            if (stream.readLine() == "          **** ALPHA SET **** ") {
                if (programm == FIREFLY) {
                    for (unsigned short i = 0; i < 5; i++)
                        stream.readLine();
                }
                return true;
            } else if (stream.readLine() == "          **** BETA SET ****") {
                return true;
            } else {
                return false;
            }
        }
    }
}

// //Лямбда парсер
// auto parsePopulationMatrix = [&](auto& targetContainer) {
//     unsigned short o = 0, a = 0, s = 0;
//     do {
//         const unsigned short block_start = o;
//         for (unsigned short j = 0; j < basis_functions + 1; ++j) {
//             QString line = stream.readLine().trimmed();
//             if (line.isEmpty()) line = stream.readLine().trimmed();
//             if (j > 0) { //Пропуск строки заголовка колонок
//                 QRegularExpressionMatch m = re.match(line);
//                 if (m.hasMatch()) {
//                     a = m.capturedView(2).toUShort() - 1; //Извлечение счетчика атома
//                     auto shell_type = m.capturedView(4);
//                     if (s != 0 && shell_type.last() == u's') {
//                         s = 0;
//                     } else if (s != 1 && shell_type[1] == u'p') {
//                         s = 1;
//                     } else if (s != 2 && shell_type[1] == u'd') {
//                         s = 2;
//                     } else if (s != 3 && shell_type[1] == u'f') {
//                         s = 3;
//                     } else if (s != 4 && shell_type[1] == u'g') {
//                         s = 4;
//                     } else if (s != 5 && shell_type[1] == u'h') {
//                         s = 5;
//                     }
//                     o = block_start;
//                     for (int i = 5; i <= m.lastCapturedIndex(); ++i) {
//                         targetContainer[o][a][s].push_back(m.capturedView(i).toDouble());
//                         S(o, j) = m.capturedView(i).toDouble();
//                         ++o;
//                     }
//                 }
//             }
//         }
//     } while (o < orbital_count);
// };
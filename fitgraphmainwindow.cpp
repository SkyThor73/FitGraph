#include "fitgraphmainwindow.h"
#include "ui_fitgraphmainwindow.h"

FitGraphMainWindow::FitGraphMainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::FitGraphMainWindow)
{
    CheckPCS();
    ui->setupUi(this);
    output = new FOutput;
    pcs = new FDataPCS;
    calculate = new FCalculation;
    exp = new FExperimental;
    adjust = new FAdjusting;
    ui->widgetPlot->xAxis->setLabel("E, eV");                   //Подпись оси x
    ui->widgetPlot->yAxis->setLabel("Intensity, a.u.");         //Подпись оси y
    ui->widgetPlot->yAxis->setTicks(false);

    shift_value = new QLabel(this);                         //Строка сдвига
    ui->statusBar->addPermanentWidget(shift_value);         //Строка состояния
    ui->widgetPlot->setInteractions(QCP::iSelectPlottables);    //Возможность выбора точек графика

    connect(ui->horizontalSlider_Bandwidth, SIGNAL(valueChanged(int)), this, SLOT(setDigitalNumber(int)));
    connect(this, SIGNAL(ChangeInitial_value(QString)), ui->lineEdit_Min, SLOT(setText(QString)));
    connect(this, SIGNAL(ChangeEnd_value(QString)), ui->lineEdit_Max, SLOT(setText(QString)));
    connect(this, SIGNAL(ChangeSteps(QString)), ui->lineEdit_Number_of_point, SLOT(setText(QString)));
    connect(this, SIGNAL(ChangeShift(QString)), ui->lineEdit_Shift, SLOT(setText(QString)));

    connect(ui->widgetPlot, SIGNAL(plottableClick(QCPAbstractPlottable*,int,QMouseEvent*)), this, SLOT(graphClicked(QCPAbstractPlottable*,int)));
    connect(ui->widgetPlot, SIGNAL(selectionChangedByUser()), this, SLOT(graphSelectionChanged()));
}

FitGraphMainWindow::~FitGraphMainWindow()
{
    delete shift_value;
    delete adjust;
    delete exp;
    delete calculate;
    delete pcs;
    delete output;
    delete ui;
}

void FitGraphMainWindow::CheckPCS()
{
    if (!QFile(FDataPCS::pcs_path).exists()) {
        int n = QMessageBox::critical(this, "System failure", "File PhotoionizationCrossSection.dat does not exist!", QMessageBox::Ok);
        if (n == QMessageBox::Ok) {
            qCritical("File '%s' does not exist!", qUtf8Printable(FDataPCS::pcs_path));
            exit(EXIT_FAILURE);
        }
    }
}

void FitGraphMainWindow::PlotCalculationGraph()
{
    ui->widgetPlot->graph(1)->data()->clear();
    ui->widgetPlot->xAxis->setRange(calculate->getInitial_value(), calculate->getEnd_value());
    if (ui->pushButton_Adjustment->isChecked()) {
        if (exp->getMax_value_exp_y() > calculate->getMax_value_y()) {
            ui->widgetPlot->yAxis->setRange(0, exp->getMax_value_exp_y());
        } else {
            ui->widgetPlot->yAxis->setRange(0, calculate->getMax_value_y());
        }
    } else {
        ui->widgetPlot->yAxis->setRange(0, calculate->getMax_value_y());
    }
    ui->widgetPlot->graph(1)->addData(calculate->getAbscissaX(), calculate->getOrdinateY());
    ui->widgetPlot->replot();
}

void FitGraphMainWindow::PlotSingleDistribution()
{
    ui->widgetPlot->graph(2)->data()->clear();
    if (selected_orbital_number != orbitals && (selected_orbital_number - 1) > 0 && (energy[selected_orbital_number - 1] == energy[selected_orbital_number - 2] || energy[selected_orbital_number - 1] == energy[selected_orbital_number])) {
        QList<double> double_gaussian = calculate->getGaussian()[selected_orbital_number - 1];
        for (auto&var : double_gaussian) {
            var += var;
        }
        ui->widgetPlot->graph(2)->addData(calculate->getAbscissaX(),  double_gaussian);
    } else {
        ui->widgetPlot->graph(2)->addData(calculate->getAbscissaX(),  calculate->getGaussian()[selected_orbital_number - 1]);   //Построение одиночного распределения выбранной МО
    }
    ui->widgetPlot->replot();
}

void FitGraphMainWindow::PlotExperimentalGraph()
{
    ui->widgetPlot->graph(0)->data()->clear();
    if (exp->getMax_value_exp_y() > calculate->getMax_value_y() || exp->getMax_value_exp_y() == calculate->getMax_value_y())
        ui->widgetPlot->yAxis->setRange(0, exp->getMax_value_exp_y());
    ui->widgetPlot->graph(0)->addData(exp->getExpAbscissaX(), exp->getExpOrdinateY());
    ui->widgetPlot->replot();
}

void FitGraphMainWindow::CreateContributionsBars()
{
    QString spd = "";                                       //Строка для имени колонки

    std::map<double, unsigned short> map;   //Класс для определения максимального количества повторяющихся элементов
    for (auto& r : energy)
        map[r]++;
    unsigned short max = 0; //Максимальное количество elementBars
    for (auto& r : map)
        if (r.second > max)
            max = r.second;
    elementBars.clear();
    elementBars.resize(max);
    for (unsigned short i = 0; i < max; ++i) {
        elementBars[i].resize(output->getUnique_atoms());
        for (unsigned short a = 0; a < output->getUnique_atoms(); a++) {
            elementBars[i][a].resize(output->getKSpdfghi());
            for (unsigned short s = 0; s < output->getKSpdfghi(); s++) {
                for (unsigned short n = 0; n < 8; n++) {
                    elementBars[i][a][s].push_back(new QCPBars(ui->widgetPlot->xAxis, ui->widgetPlot->yAxis));
                    //Распологаем колонки друг над другом
                    if (n > 0) {
                        elementBars[i][a][s][n]->moveAbove(elementBars[i][a][s][n - 1]);
                    } else if (s > 0) {
                        elementBars[i][a][s][n]->moveAbove(elementBars[i][a][s - 1][6]);
                    } else if (a > 0) {
                        elementBars[i][a][s][n]->moveAbove(elementBars[i][a - 1][output->getKSpdfghi() - 1][6]);
                    }
                    //Добавляем данные для колонок если есть для них сечения
                    if (pcs->getPhotoionizationCrossSection()[0][a][s][n] != 0 || pcs->getPhotoionizationCrossSection()[1][a][s][n] != 0 ||
                        pcs->getPhotoionizationCrossSection()[2][a][s][n] != 0 || pcs->getPhotoionizationCrossSection()[3][a][s][n] != 0) {
                        switch (s) {
                        case 0:
                            spd = "S";
                            break;
                        case 1:
                            spd = "P";
                            break;
                        case 2:
                            spd = "D";
                            break;
                        case 3:
                            spd = "F";
                            break;
                        default:
                            break;
                        }
                        elementBars[i][a][s][n]->setName(element[a]->getElementName() + '\t' + QString::number(n + s + 1) + spd + '\t' + QString::number(i));  //Имя вида: АТОМ НОМЕРподуровень
                        elementBars[i][a][s][n]->setAntialiased(false);//Выключение сглаживания
                        elementBars[i][a][s][n]->setWidthType(QCPBars::wtAbsolute);
                        elementBars[i][a][s][n]->setWidth(8);
                        elementBars[i][a][s][n]->setStackingGap(1);    //Щель между колонками в пикселях
                        //Окрас колонок
                        if (element[a]->getElementName() == "H") {
                            elementBars[i][a][s][n]->setPen(QPen(QColor(0, 0, 0)));
                        } else {
                            elementBars[i][a][s][n]->setPen(QPen(element[a]->getElementColor()));
                            elementBars[i][a][s][n]->setBrush(element[a]->getElementColor().darker(196 - 49 * s - 7 * n));
                        }
                    }
                }
            }
        }
    }
    //Располагаем все вырожденные орбитали друг над другом
    for (unsigned short i = 1; i < max; ++i) {
        elementBars[i][0][0][0]->moveAbove(elementBars[i - 1][output->getUnique_atoms() - 1][output->getKSpdfghi() - 1][6]);
        for (unsigned short a = 0; a < output->getUnique_atoms(); a++) {
            for (unsigned short s = 0; s < output->getKSpdfghi(); s++) {
                for (unsigned short n = 0; n < 8; n++) {
                    //Распологаем колонки друг над другом
                    if (n > 0) {
                        elementBars[i][a][s][n]->moveAbove(elementBars[i][a][s][n - 1]);
                    } else if (s > 0) {
                        elementBars[i][a][s][n]->moveAbove(elementBars[i][a][s - 1][6]);
                    } else if (a > 0) {
                        elementBars[i][a][s][n]->moveAbove(elementBars[i][a - 1][output->getKSpdfghi() - 1][6]);
                    }
                }
            }
        }
    }
}

void FitGraphMainWindow::PlotContributionsBars()
{
    QList<QList<double>> temp_ticks_energy_vec;//Данные с энергиями для колонок
    QList<QList<double>> temp_values_vec;           //Данные со значениями для колонок
    temp_ticks_energy_vec.resize(elementBars.size());
    temp_values_vec.resize(elementBars.size());
    for (unsigned short a = 0; a < output->getUnique_atoms(); a++) {
        for (unsigned short s = 0; s < output->getKSpdfghi(); s++) {
            for (unsigned short n = 0; n < 8; n++) {
                for (unsigned short i = 0; i < elementBars.size(); ++i) {
                    temp_ticks_energy_vec[i].clear();
                    temp_values_vec[i].clear();
                }
                unsigned short i = 0;
                for (unsigned short o = 0; o < orbitals; o++) {
                    if (o > 0 && energy[o] == energy[o - 1]) {
                        ++i;
                        if (calculate->getSpdRIntensity()[o][a][s][n] != 0) {
                            temp_ticks_energy_vec[i].push_back(energy[o]);
                            temp_values_vec[i].push_back(calculate->getSpdRIntensity()[o][a][s][n]);
                        }
                    } else {
                        i = 0;
                        if (calculate->getSpdRIntensity()[o][a][s][n] != 0) {
                            temp_ticks_energy_vec[i].push_back(energy[o]);
                            temp_values_vec[i].push_back(calculate->getSpdRIntensity()[o][a][s][n]);
                        }
                    }
                }
                for (unsigned short i = 0; i < elementBars.size(); ++i) {
                    elementBars[i][a][s][n]->data()->clear();
                    elementBars[i][a][s][n]->setData(temp_ticks_energy_vec[i], temp_values_vec[i]);
                }
            }
        }
    }
    ui->widgetPlot->replot();
}

void FitGraphMainWindow::on_pushButton_Open_clicked()
{
    setPath_file(QFileDialog::getOpenFileName(this, "Choose file", "", "*.out *.log"));//Попытка присвоить путь к файлу
    if (getPath_file() != "") {                             //Если путь найден, программа открывает файл и работает с ним
        ui->widgetPlot->clearGraphs();                          //Очищаем все данные обо всех графиках
        ui->widgetPlot->addGraph();
        QPen graphPen0;
        graphPen0.setWidthF(2);
        graphPen0.setColor(QColor(255, 0, 0));
        ui->widgetPlot->graph(0)->setPen(graphPen0);
        ui->widgetPlot->graph(0)->setName("Experimental");
        ui->widgetPlot->addGraph();
        QPen graphPen1;
        graphPen1.setWidthF(2);
        graphPen1.setColor(QColor(0, 0, 255));
        ui->widgetPlot->graph(1)->setPen(graphPen1);
        ui->widgetPlot->graph(1)->setName("Calculated");
        ui->widgetPlot->addGraph();
        ui->widgetPlot->graph(2)->setPen(QColor(0, 0, 0));
        ui->widgetPlot->graph(2)->setSelectable(QCP::SelectionType::stNone);
        ui->widgetPlot->graph(2)->setName("MO");
        //Сброс нажатых кнопок
        ui->pushButton_Adjustment->setChecked(false);
        ui->pushButton_Contributions->setChecked(false);
        //Установка радиобаттонов по умолчанию
        ui->radioButton_He_I->setChecked(true);
        ui->radioButton_Simple->setChecked(true);
        ui->radioButton_Alpha->setChecked(true);

        output->ReadData(getPath_file());                   //Считывание всех необходимых данных из выходного файла

        calculate->setInitial_value(0);                     //Установка начального положения графика
        emit ChangeInitial_value(QString::number(calculate->getInitial_value()));
        calculate->setResolution_ratio(1.0 * (ui->horizontalSlider_Bandwidth->value()) / 100);//Установка множителя интенсивности
        calculate->setSteps(1000);                          //Установка количества точек
        emit ChangeSteps(QString::number(calculate->getSteps()));
        calculate->setSource(FCalculation::HE_I);           //Установка сечения
        calculate->setUpper_bound(20);                      //Установка верхнего предела на ионизацию
        calculate->setEnd_value(20);                        //Установка конечного положения графика
        emit ChangeEnd_value(QString::number(calculate->getEnd_value()));
        //Создание объектов класса элементов
        element.clear();
        for (unsigned short a = 0; a < output->getUnique_atoms(); a++) {
            element.push_back(new FElement);
        }
        pcs->MiningPCS(*output, element);                   //Сбор сечений фотоионизаций нужных атомов и заполнение данными элементов

        calculate->Calculation(*output, element);           //Вычисление вкладов

        orbitals = output->getOrbital_count();              //Присваивание переменной orbitals количество орбиталей альфа набора по умолчанию
        energy = output->getEnergy();                       //Присваивание переменной energy альфа набор энергий
        spdVectors = calculate->getSpdSimple();             //Присвоение вектору для построения набор простых вкладов по умолчанию

        calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);//Вычисление интенсивностей и графика

        std::ranges::sort(element, std::greater(), &FElement::getElementNumber);
        CreateContributionsBars();                          //Создание объектов вкладов
        //ui->widgetPlot->setToolTip("name");                   //TODO:показ брутто формулы Stoichiometry RUN TITLE
        //Включение/выключение кнопок
        ui->radioButton_He_I->setEnabled(true);
        ui->radioButton_He_II->setEnabled(true);
        ui->radioButton_Al->setEnabled(true);
        ui->radioButton_Mg->setEnabled(true);
        ui->horizontalSlider_Bandwidth->setEnabled(true);
        ui->lineEdit_Max->setEnabled(true);
        ui->lineEdit_Min->setEnabled(true);
        ui->lineEdit_Number_of_point->setEnabled(true);

        ui->pushButton_Save->setEnabled(true);
        ui->pushButton_Add->setEnabled(true);
        ui->pushButton_Contributions->setEnabled(true);
        //ui->pushButton_Contributions->setCheckable(true);
        ui->pushButton_Adjustment->setEnabled(false);
        ui->radioButton_Simple->setEnabled(true);
        if (output->getUhf()) {
            ui->radioButton_Alpha->setEnabled(true);
            ui->radioButton_Beta->setEnabled(true);
        } else {
            ui->radioButton_Alpha->setEnabled(false);
            ui->radioButton_Beta->setEnabled(false);
        }
        if (output->getLowdin()) {
            ui->radioButton_Lowdin->setEnabled(true);
        } else {
            ui->radioButton_Lowdin->setEnabled(false);
        }
        if (output->getMulliken()) {
            ui->radioButton_Mulliken->setEnabled(true);
        } else {
            ui->radioButton_Mulliken->setEnabled(false);
        }
        exp->setMax_value_exp_y(0);
        PlotCalculationGraph();
    }
}


void FitGraphMainWindow::on_pushButton_Add_clicked()
{
    exp->setPath_exp_file(QFileDialog::getOpenFileName(this, "Choose file", "", "*.dat"));//Путь файла с экспериментальным спектром
    if (exp->getPath_exp_file() != "") {
        exp->ReadXY(*calculate);  //Чтение и сбор данных

        calculate->setInitial_value(exp->getExpAbscissaX().front());//Установка начального значения в соответствии с начальным значением экспериментального спектра
        emit ChangeInitial_value(QString::number(calculate->getInitial_value()));
        calculate->setEnd_value(exp->getExpAbscissaX().back()); //Установка конечного значения в соответствии с конечным значением экспериментального спектра
        emit ChangeEnd_value(QString::number(calculate->getEnd_value()));
        //calculate->setUpper_bound(exp->getExpAbscissaX().back());
        calculate->setSteps(exp->getExpAbscissaX().size());     //Установка точек графика в соответствии с точками экспериментального спектра
        emit ChangeSteps(QString::number(calculate->getSteps()));
        calculate->CalculationXY(orbitals, energy);             //Пересчёт области спектра
        //При изменении максимального значения расчётного спектра
        //double multiplier = 0;
        //double temp_sum_I_calc = 0;
        // for (unsigned short i = 4; i > 0; --i)
        //      temp_sum_I_calc += calculate->getIntegral_area()[i];
        //temp_sum_I_calc += calculate->getIntegral_area()[0];
        //multiplier = temp_sum_I_calc / 637.4;
        double multiplier = calculate->getMax_value_y() / exp->getMax_value_exp_y();
        QList<double> experimentalY = exp->getExpOrdinateY();
        for (int i = 0; i < exp->getExpOrdinateY().size(); ++i) {
            experimentalY[i] *= multiplier;
        }
        exp->setExpOrdinateY(experimentalY);
        exp->setMax_value_exp_y(*std::max_element(experimentalY.begin(), experimentalY.end()));

        PlotCalculationGraph();
        if (ui->pushButton_Contributions->isChecked())
            PlotContributionsBars();
        PlotExperimentalGraph();                                //Построение экспериментального графика
        ui->pushButton_Adjustment->setEnabled(true);                //Разблокирование кнопки подгона
    }
}


void FitGraphMainWindow::on_pushButton_Save_clicked()
{
    QString strFilter;
    QString file_name = getPath_file();
    file_name.chop(4);
    QString str = QFileDialog::getSaveFileName(this, tr("Save graph"), file_name, "*.txt ;; *.dat" , &strFilter);
    QFile fileOut(str);
    if (fileOut.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream writeStream(&fileOut);
        for (unsigned short x = 0; x < calculate->getSteps(); x++) {
            writeStream << calculate->getOrdinateY()[x] << '\t' << calculate->getAbscissaX()[x] << '\n';
        }
        fileOut.close();
    }
}


void FitGraphMainWindow::on_radioButton_He_I_clicked()
{
    calculate->setSource(FCalculation::HE_I);
    calculate->setInitial_value(0);
    emit ChangeInitial_value(QString::number(calculate->getInitial_value()));
    calculate->setResolution_ratio(1.0 * (ui->horizontalSlider_Bandwidth->value()) / 100);
    calculate->setSteps(1000);
    emit ChangeSteps(QString::number(calculate->getSteps()));
    calculate->setUpper_bound(20);
    calculate->setEnd_value(20);
    emit ChangeEnd_value(QString::number(calculate->getEnd_value()));
    calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
}


void FitGraphMainWindow::on_radioButton_He_II_clicked()
{
    calculate->setSource(FCalculation::HE_II);
    calculate->setInitial_value(0);
    emit ChangeInitial_value(QString::number(calculate->getInitial_value()));
    calculate->setResolution_ratio(1.0 * (ui->horizontalSlider_Bandwidth->value()) / 100);
    calculate->setSteps(1000);
    emit ChangeSteps(QString::number(calculate->getSteps()));
    calculate->setUpper_bound(40);
    calculate->setEnd_value(40);
    emit ChangeEnd_value(QString::number(calculate->getEnd_value()));
    calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
}


void FitGraphMainWindow::on_radioButton_Mg_clicked()
{
    calculate->setSource(FCalculation::MG);
    calculate->setInitial_value(0);
    emit ChangeInitial_value(QString::number(calculate->getInitial_value()));
    calculate->setResolution_ratio(1.0 * (ui->horizontalSlider_Bandwidth->value()) / 100);
    calculate->setSteps(2000);
    emit ChangeSteps(QString::number(calculate->getSteps()));
    calculate->setUpper_bound(1200);
    calculate->setEnd_value(1200);
    emit ChangeEnd_value(QString::number(calculate->getEnd_value()));
    calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
}


void FitGraphMainWindow::on_radioButton_Al_clicked()
{
    calculate->setSource(FCalculation::AL);
    calculate->setInitial_value(0);
    emit ChangeInitial_value(QString::number(calculate->getInitial_value()));
    calculate->setResolution_ratio(1.0 * (ui->horizontalSlider_Bandwidth->value()) / 100);
    calculate->setSteps(2000);
    emit ChangeSteps(QString::number(calculate->getSteps()));
    calculate->setUpper_bound(1400);
    calculate->setEnd_value(1400);
    emit ChangeEnd_value(QString::number(calculate->getEnd_value()));
    calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
}


void FitGraphMainWindow::on_radioButton_Simple_clicked()
{
    if (ui->radioButton_Alpha->isChecked()) {
        spdVectors = calculate->getSpdSimple();
    } else if (ui->radioButton_Beta->isChecked()) {
        spdVectors = calculate->getSpdSimpleBeta();
    }

    calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
}


void FitGraphMainWindow::on_radioButton_Mulliken_clicked()
{
    if (ui->radioButton_Alpha->isChecked()) {
        spdVectors = calculate->getSpdMulliken();
    } else if (ui->radioButton_Beta->isChecked()) {
        spdVectors = calculate->getSpdMullikenBeta();
    }

    calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
}


void FitGraphMainWindow::on_radioButton_Lowdin_clicked()
{
    if (ui->radioButton_Alpha->isChecked()) {
        spdVectors = calculate->getSpdLowdin();
    } else if (ui->radioButton_Beta->isChecked()) {
        spdVectors = calculate->getSpdLowdinBeta();
    }
    spdVectors = calculate->getSpdLowdin();

    calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
}


void FitGraphMainWindow::on_radioButton_Alpha_clicked()
{
    orbitals = output->getOrbital_count();                  //Присвоение переменной orbitals количество орбиталей альфа набора
    energy = output->getEnergy();
    if (ui->radioButton_Simple->isChecked()) {
        spdVectors = calculate->getSpdSimple();
    } else if (ui->radioButton_Mulliken->isChecked()) {
        spdVectors = calculate->getSpdMulliken();
    } else if (ui->radioButton_Lowdin->isChecked()) {
        spdVectors = calculate->getSpdLowdin();
    }

    calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);

    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
    PlotCalculationGraph();
}


void FitGraphMainWindow::on_radioButton_Beta_clicked()
{
    orbitals = output->getOrbital_count_beta();             //Присвоение переменной orbitals количество орбиталей бета набора
    energy = output->getEnergyBeta();
    if (ui->radioButton_Simple->isChecked()) {
        spdVectors = calculate->getSpdSimpleBeta();
    } else if (ui->radioButton_Mulliken->isChecked()) {
        spdVectors = calculate->getSpdMullikenBeta();
    } else if (ui->radioButton_Lowdin->isChecked()) {
        spdVectors = calculate->getSpdLowdinBeta();
    }

    calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);

    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
    PlotCalculationGraph();
}


void FitGraphMainWindow::on_pushButton_Contributions_toggled(bool checked)
{
    if (checked) {
        PlotContributionsBars();
        PlotCalculationGraph();
    } else {
        for (unsigned short a = 0; a < output->getUnique_atoms(); a++) {
            for (unsigned short s = 0; s < output->getKSpdfghi(); s++) {
                for (unsigned short n = 0; n < 8; n++) {
                    for (int i = 0; i < elementBars.size(); ++i) {
                        elementBars[i][a][s][n]->data()->clear();
                    }
                }
            }
        }
        PlotCalculationGraph();
    }
    if (ui->pushButton_Adjustment->isEnabled() && ui->pushButton_Adjustment->isChecked())
        PlotExperimentalGraph();
}


void FitGraphMainWindow::on_pushButton_Adjustment_toggled(bool checked)
{
    if (checked) {
        calculate->setInitial_value(exp->getExpAbscissaX().front());//Установка начального значения в соответствии с начальным значением экспериментального спектра
        emit ChangeInitial_value(QString::number(calculate->getInitial_value()));
        calculate->setEnd_value(exp->getExpAbscissaX().back());//Установка конечного значения в соответствии с конечным значением экспериментального спектра
        emit ChangeEnd_value(QString::number(calculate->getEnd_value()));
        calculate->setSteps(exp->getExpAbscissaX().size()); //Установка точек графика в соответствии с точками экспериментального спектра
        emit ChangeSteps(QString::number(calculate->getSteps()));
        //Пересчёт области спектра
        QList<double> temp_resolution_vector;
        temp_resolution_vector.resize(orbitals);
        // if (calculate->getSource() == FCalculation::HE_I) {
        //     for (int o = 0; o < orbitals; ++o) {
        //         temp_resolution_vector[o] = 0.1741 * pow(M_E, 0.1135 * energy[o]); // 0.2305 * energy[o] - 2.5004
        //     }
        // } else {
        //     for (auto& o : temp_resolution_vector) {
        //         o = calculate->getResolution_ratio();
        //     }
        // }
        for (auto& o : temp_resolution_vector) {
            o = calculate->getResolution_ratio();
        }
        calculate->setResolution_vector(temp_resolution_vector);
        calculate->CalculationXY(orbitals, energy, calculate->getResolution_vector());
        PlotCalculationGraph();
        if (ui->pushButton_Contributions->isChecked())
            PlotContributionsBars();
        PlotExperimentalGraph();
        //Установка кнопок
        ui->lineEdit_Max->setEnabled(false);
        ui->lineEdit_Min->setEnabled(false);
        ui->lineEdit_Number_of_point->setEnabled(false);
        ui->radioButton_Al->setEnabled(false);
        ui->radioButton_He_I->setEnabled(false);
        ui->radioButton_He_II->setEnabled(false);
        ui->radioButton_Mg->setEnabled(false);
        ui->radioButton_Mulliken->setEnabled(false);
        ui->radioButton_Simple->setEnabled(false);
        ui->radioButton_Lowdin->setEnabled(false);
        ui->radioButton_Alpha->setEnabled(false);
        ui->radioButton_Beta->setEnabled(false);
        ui->pushButton_Shift->setEnabled(true);
        ui->lineEdit_Shift->setEnabled(true);
        ui->pushButton_Vise->setEnabled(true);
        ui->pushButton_Fit->setEnabled(true);
    } else {
        //Возврат значений энергий МО к изначальным
        if (ui->radioButton_Alpha->isChecked()) { //TODO
            energy = output->getEnergy();
        } else if (ui->radioButton_Beta->isChecked()) {
            energy = output->getEnergyBeta();
        }
        adjust->setVisible_shift_x(0.0);
        //Пересчёт области спектра
        ui->widgetPlot->graph(0)->data()->clear();
        calculate->CalculationXY(orbitals, energy);
        PlotCalculationGraph();
        if (ui->pushButton_Contributions->isChecked())
            PlotContributionsBars();
        //Установка кнопок
        ui->pushButton_Shift->setEnabled(false);
        ui->lineEdit_Shift->setEnabled(false);
        ui->pushButton_Vise->setEnabled(false);
        ui->pushButton_Fit->setEnabled(false);
        ui->lineEdit_Max->setEnabled(true);
        ui->lineEdit_Min->setEnabled(true);
        ui->lineEdit_Number_of_point->setEnabled(true);
        ui->radioButton_Al->setEnabled(true);
        ui->radioButton_He_I->setEnabled(true);
        ui->radioButton_He_II->setEnabled(true);
        ui->radioButton_Mg->setEnabled(true);
        ui->radioButton_Simple->setEnabled(true);
        if (output->getUhf()) {
            ui->radioButton_Alpha->setEnabled(true);
            ui->radioButton_Beta->setEnabled(true);
        } else {
            ui->radioButton_Alpha->setEnabled(false);
            ui->radioButton_Beta->setEnabled(false);
        }
        if (output->getLowdin()) {
            ui->radioButton_Lowdin->setEnabled(true);
        } else {
            ui->radioButton_Lowdin->setEnabled(false);
        }
        if (output->getMulliken()) {
            ui->radioButton_Mulliken->setEnabled(true);
        } else {
            ui->radioButton_Mulliken->setEnabled(false);
        }
        shift_value->clear();
    }
}


void FitGraphMainWindow::on_pushButton_Shift_clicked()
{
    adjust->FindBetterShift(orbitals, energy, *calculate, *exp);
    adjust->ShiftEnergy(energy);
    calculate->CalculationXY(orbitals, energy, calculate->getResolution_vector());

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
    PlotExperimentalGraph();

    adjust->CalculationCompliance(*calculate, *exp);
    shift_value->setText("Δε = " + QString().setNum(round(adjust->getVisible_shift_x()*10000)/10000) + " eV" + '\t' + "Matching " + QString().setNum(adjust->getSpecific_min_square()));
}


void FitGraphMainWindow::on_pushButton_Vise_clicked()
{
    adjust->FindBetterVise(orbitals, energy, *calculate, *exp);

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
    PlotExperimentalGraph();

    adjust->CalculationCompliance(*calculate, *exp);
    shift_value->setText("Δε = " + QString().setNum(round(adjust->getVisible_shift_x()*10000)/10000) + " eV" + '\t' + "Matching " + QString().setNum(adjust->getSpecific_min_square()));
}


void FitGraphMainWindow::on_pushButton_Fit_clicked()
{
    adjust->FindBetterWidth(orbitals, energy, *calculate, *exp);
    calculate->CalculationEachWidth(*output, *pcs, orbitals, spdVectors, calculate->getResolution_vector());

    PlotCalculationGraph();
    PlotExperimentalGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();

    adjust->CalculationCompliance(*calculate, *exp);
    shift_value->setText("Δε = " + QString().setNum(round(adjust->getVisible_shift_x()*10000)/10000) + " eV" + '\t' + "Matching " + QString().setNum(adjust->getSpecific_min_square()));
}


void FitGraphMainWindow::on_lineEdit_Shift_editingFinished()
{
    energy = output->getEnergy();
    adjust->setVisible_shift_x((ui->lineEdit_Shift->text()).toDouble());
    emit ChangeShift(QString::number(adjust->getVisible_shift_x()));
    for (unsigned short o = 0; o < energy.size(); o++)
        energy[o] += adjust->getVisible_shift_x();
    calculate->CalculationXY(orbitals, energy, calculate->getResolution_vector());

    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
    PlotExperimentalGraph();

    adjust->CalculationCompliance(*calculate, *exp);
    shift_value->setText("Δε = " + QString().setNum(round(adjust->getVisible_shift_x()*10000)/10000) + " eV" + '\t' + "Matching " + QString().setNum(adjust->getSpecific_min_square()));
}


void FitGraphMainWindow::on_checkBox_Transparent_stateChanged(int)
{
    if (ui->checkBox_Transparent->isChecked()) {
        setWindowOpacity(0.7);
    } else {
        setWindowOpacity(1.0);
    }
}


void FitGraphMainWindow::on_lineEdit_Number_of_point_editingFinished()
{
    calculate->setSteps((ui->lineEdit_Number_of_point->text()).toUShort());
    emit ChangeSteps(QString::number(calculate->getSteps()));
    calculate->CalculationXY(orbitals, energy);

    PlotCalculationGraph();
}


void FitGraphMainWindow::on_lineEdit_Min_editingFinished()
{
    calculate->setInitial_value((ui->lineEdit_Min->text()).toDouble());
    emit ChangeInitial_value(QString::number(calculate->getInitial_value()));
    calculate->CalculationXY(orbitals, energy);

    PlotCalculationGraph();
}


void FitGraphMainWindow::on_lineEdit_Max_editingFinished()
{
    calculate->setEnd_value((ui->lineEdit_Max->text()).toDouble());
    emit ChangeEnd_value(QString::number(calculate->getEnd_value()));
    calculate->CalculationXY(orbitals, energy);

    PlotCalculationGraph();
}


void FitGraphMainWindow::on_horizontalSlider_Bandwidth_valueChanged(int value)
{
    if (ui->pushButton_Adjustment->isChecked()) {
        if (check_selected) {
            QList<double> temp_resolution_vector = calculate->getResolution_vector();
            if (energy[selected_orbital_number - 1] == energy[selected_orbital_number - 2]) {
                temp_resolution_vector[selected_orbital_number - 2] = 1.0 * value / 100;
            } else if(selected_orbital_number != orbitals && energy[selected_orbital_number - 1] == energy[selected_orbital_number]) {
                temp_resolution_vector[selected_orbital_number] = 1.0 * value / 100;
            }
            temp_resolution_vector[selected_orbital_number - 1] = 1.0 * value / 100;
            calculate->setResolution_vector(temp_resolution_vector);
            calculate->CalculationXY(orbitals, energy, calculate->getResolution_vector());
            calculate->CalculationEachWidth(*output, *pcs, orbitals, spdVectors, calculate->getResolution_vector());
            adjust->CalculationCompliance(*calculate, *exp);
            shift_value->setText("Δε = " + QString().setNum(round(adjust->getVisible_shift_x()*10000)/10000) + " eV" + '\t' + "Matching " + QString().setNum(adjust->getSpecific_min_square()));
        } else {
            calculate->setResolution_ratio(1.0 * value / 100);
            //calculate->CalculationXY(orbitals, energy);
            calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);


            QList<double> temp_resolution_vector;
            temp_resolution_vector.resize(orbitals);
            for (auto& o : temp_resolution_vector) {
                o = calculate->getResolution_ratio();
            }

            calculate->setResolution_vector(temp_resolution_vector);
            adjust->CalculationCompliance(*calculate, *exp);
            shift_value->setText("Δε = " + QString().setNum(round(adjust->getVisible_shift_x()*10000)/10000) + " eV" + '\t' + "Matching " + QString().setNum(adjust->getSpecific_min_square()));
        }
        PlotExperimentalGraph();
    } else {
        calculate->setResolution_ratio(1.0 * value / 100);
        calculate->CalculationIntegralArea(*output, *pcs, orbitals, energy, spdVectors);
    }
    PlotCalculationGraph();
    if (ui->pushButton_Contributions->isChecked())
        PlotContributionsBars();
    if (check_selected) {
        PlotSingleDistribution();
    }
}

void FitGraphMainWindow::setDigitalNumber(int value)
{
    ui->label_Bandwidth_Count->setNum(1.0 * value / 100);
}

void FitGraphMainWindow::graphSelectionChanged()
{
    for (int i=0; i<ui->widgetPlot->graphCount(); ++i)
    {
        QCPGraph *graph = ui->widgetPlot->graph(i);
        QCPPlottableLegendItem *item = ui->widgetPlot->legend->itemWithPlottable(graph);
        if (item->selected() || graph->selected()) {
            item->setVisible(true);
            graph->setSelection(QCPDataSelection(graph->data()->dataRange()));
        }
    }
    ui->widgetPlot->graph(2)->data()->clear();      //Очистка данных одиночного распределения, при нажатии на любую точку окна построения
    check_selected = false;
}

void FitGraphMainWindow::graphClicked(QCPAbstractPlottable *plottable, int dataIndex)
{
    if (plottable->name() != "MO") {
        double dataKey = plottable->interface1D()->dataMainKey(dataIndex);  //Значения энергии МО
        QString message = "";
        if (plottable->name() == "Calculated" || plottable->name() == "Experimental") {
            message = QString("%1 value %2 energy %3").arg(plottable->name()).arg(plottable->interface1D()->dataMainValue(dataIndex)).arg(dataKey);
        } else {
            check_selected = true;
            QString dataName = "";
            for (unsigned short o = orbitals - 1; o >= 0; --o) {
                if (dataKey == energy[o]) {
                    o -= plottable->name().last(1).toUShort();
                    selected_orbital_number = o + 1;
                    if ((o == (orbitals - 1) && (energy[o] == energy[o - 1] && output->getSymmetry()[o] == output->getSymmetry()[o - 1])) || (o == 0 && (energy[o] == energy[o + 1] && output->getSymmetry()[o] == output->getSymmetry()[o + 1]))) {
                        dataName = plottable->name().first(5) + '\t' + output->getSymmetry()[o] + '*';
                    } else if (o != 0 && ((energy[o] == energy[o - 1] && output->getSymmetry()[o] == output->getSymmetry()[o - 1]) || (o != (orbitals - 1) && (energy[o] == energy[o + 1] && output->getSymmetry()[o] == output->getSymmetry()[o + 1])))) {
                        dataName = plottable->name().first(5) + '\t' + output->getSymmetry()[o] + '*';
                    } else {
                        dataName = plottable->name().first(5) + '\t' + output->getSymmetry()[o];
                    }
                    break;
                }
            }
            double dataValue = round(plottable->interface1D()->dataMainValue(dataIndex) / calculate->getRIntensity()[selected_orbital_number - 1] * 100);
            if (ui->pushButton_Adjustment->isChecked()) {
                message = QString("%1 %2 contribution %3 energy %4 width %5 area %6").arg(selected_orbital_number).arg(dataName).arg(dataValue).arg(dataKey).arg(calculate->getResolution_vector()[selected_orbital_number - 1]).arg(calculate->getIntegral_area()[selected_orbital_number - 1]);
            } else {
                message = QString("%1 %2 contribution %3 energy %4").arg(selected_orbital_number).arg(dataName).arg(dataValue).arg(dataKey);
            }

            PlotSingleDistribution();
        }
        ui->statusBar->showMessage(message);                    //Отображение в статус выбранного значения графика
    }
}

const QString &FitGraphMainWindow::getPath_file() const
{
    return path_file;
}

void FitGraphMainWindow::setPath_file(const QString &newPath_file)
{
    path_file = newPath_file;
}


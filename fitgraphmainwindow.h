#ifndef FITGRAPHMAINWINDOW_H
#define FITGRAPHMAINWINDOW_H

#include <QMainWindow>
#include "fadjusting.h"
#include "fcalculation.h"
#include "fdatapcs.h"
#include "fexperimental.h"
#include "foutput.h"
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class FitGraphMainWindow;
}
QT_END_NAMESPACE

class FitGraphMainWindow : public QMainWindow
{
    Q_OBJECT

public:
    FitGraphMainWindow(QWidget *parent = nullptr);
    ~FitGraphMainWindow();

    void CheckPCS();                        //Проверка существования файла PhotoionizationCrossSection.dat в папке с программой
    void PlotCalculationGraph();            //Выводит график спектра в окно программы
    void PlotSingleDistribution();
    void PlotExperimentalGraph();           //Выводит график экспериментального спектра в окно программы
    void CreateContributionsBars();         //Создаёт колонки вкладов
    void PlotContributionsBars();           //Выводит столбики вкладов в окно программы

    const QString &getPath_file() const;
    void setPath_file(const QString &newPath_file);

private slots:
    //Открывает файл и строит график по полученным данным из него
    void on_pushButton_Open_clicked();
    //Добавляет экспериментальный спектр
    void on_pushButton_Add_clicked();
    //Сохраняет график в виде файла с координатами точек графика
    void on_pushButton_Save_clicked();
    //Строит график по сечениям He I
    void on_radioButton_He_I_clicked();
    //Строит график по сечениям He II
    void on_radioButton_He_II_clicked();
    //Строит график по сечениям Mg
    void on_radioButton_Mg_clicked();
    //Строит график по сечениям Al
    void on_radioButton_Al_clicked();
    //Строит график по простым вкладам
    void on_radioButton_Simple_clicked();
    //Строит график по маликеновским вкладам
    void on_radioButton_Mulliken_clicked();
    //Строит график по вкладам Лёвдина
    void on_radioButton_Lowdin_clicked();
    //Строит график альфа набора орбиталей
    void on_radioButton_Alpha_clicked();
    //Строит график бета набора орбиталей
    void on_radioButton_Beta_clicked();
    //Включает/отключает столбики с вкладами
    void on_pushButton_Contributions_toggled(bool checked);
    //Включает/отключает режим наложения спектров
    void on_pushButton_Adjustment_toggled(bool checked);
    //Сдвигает расчётный спектр под эксперимент
    void on_pushButton_Shift_clicked();
    //Приравнивает интегральные площади спектров
    void on_pushButton_Vise_clicked();
    //Подгоняет ширину полос теоретического спектра под эксперимент
    void on_pushButton_Fit_clicked();
    //Устанавливает значение сдвига расчётного спектра
    void on_lineEdit_Shift_editingFinished();
    //Устанавлвает прозрачность окна программы (70%)
    void on_checkBox_Transparent_stateChanged(int);
    //Устанавливает значение точек графика
    void on_lineEdit_Number_of_point_editingFinished();
    //Устанавливает значение начала графика
    void on_lineEdit_Min_editingFinished();
    //Устанавливает значение конца графика
    void on_lineEdit_Max_editingFinished();
    //Меняет резкость пиков
    void on_horizontalSlider_Bandwidth_valueChanged(int value);
    //Устанавливает значение множителя резкости в окошко под ползунком
    void setDigitalNumber(int value);
    //Изменяет выбор графика
    void graphSelectionChanged();
    //Выбирает график
    void graphClicked(QCPAbstractPlottable *plottable, int dataIndex);

signals:
    //Изменение числа точек графика
    void ChangeSteps(QString);
    //Изменение начального значения графика
    void ChangeInitial_value(QString);
    //Изменение конечного значения графика
    void ChangeEnd_value(QString);
    //Изменение значения сдвига
    void ChangeShift(QString);

private:
    Ui::FitGraphMainWindow *ui;
    FOutput *output;
    FDataPCS *pcs;
    FCalculation *calculate;
    FExperimental *exp;
    FAdjusting *adjust;
    QLabel *shift_value;
    QList<QList<QList<QList<QCPBars*> > > > elementBars;
    QList<FElement*> element;

    QString path_file = "";
    unsigned short orbitals = 0;
    unsigned short selected_orbital_number = 0;
    bool check_selected = false;
    QList<double> energy;
    QList<QList<QList<QList<double> > > > spdVectors;
};
#endif // FITGRAPHMAINWINDOW_H

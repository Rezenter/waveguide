#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <complex>
#include <QDebug>
#include <QFile>
#include <QTextStream>
#include <QLineSeries>
#include <QChartView>
#include <QValueAxis>
#include <QCoreApplication>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    void calc();
    std::complex<double> i = -1;
    std::complex<double> gp;
    std::complex<double> kr;
    std::complex<double> gs;
    std::complex<double> ga;
    double na = 1.0;
    double nr = 2.4;
    double np = 1.5;
    double ns = 1.45;
    double hs = 100;
    double hr1 = 200;
    double hr2 = 200;
    double l;//lambda
    double b;//beta
    double eps;
    double omega;
    double c0 = 299792458;
    QFile *file;
    QTextStream *stream;
    std::complex<double> coeffs[800][10];
    std::complex<double> convert(QString);
    void deb(std::complex<double>);
    QtCharts::QLineSeries *zero;
    QtCharts::QLineSeries *Hxi;
    QtCharts::QLineSeries *Eyi;
    QtCharts::QLineSeries *Hxr;
    QtCharts::QLineSeries *Eyr;
    QtCharts::QLineSeries *Ezr;
    QtCharts::QLineSeries *Ezi;
    QtCharts::QLineSeries *Hxa;
    QtCharts::QLineSeries *Eya;
    QtCharts::QLineSeries *Hxm;
    QtCharts::QLineSeries *Eym;
    QtCharts::QLineSeries *Eza;
    QtCharts::QLineSeries *Ezm;
    QtCharts::QChart *chart;
    QtCharts::QValueAxis *axisX;
    QtCharts::QValueAxis *axisY;
    QtCharts::QChartView *chartView;
    void resizeEvent(QResizeEvent*);
    QFile *outFile;
    QTextStream *outStream;

private slots:
    void save();
};

#endif // MAINWINDOW_H

#include "mainwindow.h"
#include "ui_mainwindow.h"
#define k (2*M_PI/l)

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    i = sqrt(i);
    ui->setupUi(this);
    QObject::connect(ui->saveButton, &QPushButton::pressed, this, &MainWindow::save);
    QObject::connect(ui->buttonGroup, static_cast<void(QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked), this, &MainWindow::calc);
    QObject::connect(ui->switchGroup, static_cast<void(QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked), this, &MainWindow::calc);
    QObject::connect(ui->pointSpinBox, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged), this, &MainWindow::calc);
    QObject::connect(ui->lSlider, static_cast<void(QSlider::*)(int)>(&QSlider::valueChanged), this, &MainWindow::calc);
    file = new QFile("precalc.txt");
    file->open(QIODevice::ReadOnly);
    stream = new QTextStream(file);
    int c  = 0;
    setWindowTitle("waveguide");
    ui->lineEdit->setText(QCoreApplication::applicationDirPath() + "/waweguide.txt");
    while(!stream->atEnd()) {
        QString line = stream->readLine();
        if(line.length() > 0){
            QStringList args = line.split("\t");
            coeffs[c][0] = convert(args.at(0));
            coeffs[c][1] = convert(args.at(1));
            coeffs[c][2] = convert(args.at(2));
            coeffs[c][3] = convert(args.at(3));
            coeffs[c][4] = convert(args.at(4));
            coeffs[c][5] = convert(args.at(5));
            coeffs[c][6] = convert(args.at(6));
            coeffs[c][7] = convert(args.at(7));
            coeffs[c][8] = convert(args.at(8));
            coeffs[c][9] = convert(args.at(9));
        }
        c++;
    }
    file->close();

    ui->naLabel->setText(QString::number(na));
    ui->nrLabel->setText(QString::number(nr));
    ui->nsLabel->setText(QString::number(ns));
    ui->nprLabel->setText(QString::number(np));
    ui->hr1Label->setText(QString::number(hr1));
    ui->hsLabel->setText(QString::number(hs));
    ui->hr2Label->setText(QString::number(hr2));
    hs = hs*pow(10, -9);
    hr1 = hr1*pow(10, -9);
    hr2 = hr2*pow(10, -9);
    eps = 8.85*pow(10, -12);

    chart = new QtCharts::QChart();
    axisX = new QtCharts::QValueAxis;
    axisX->setMinorGridLineVisible(true);
    axisY = new QtCharts::QValueAxis();
    axisY->setTickCount(15);
    axisY->setTitleText("y[m]");
    axisX->setTickCount(10);

    ui->chartWidget->resize(450,400);
    chartView = new QtCharts::QChartView(chart, ui->chartWidget);
    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->show();
    calc();
}

MainWindow::~MainWindow()
{
    QtCharts::QAbstractSeries *delTmp;
    foreach(delTmp, chart->series()){
        delete &delTmp;
    }
    delete chart;
    delete chartView;
    delete ui;
}

std::complex<double> MainWindow::convert(QString s){
    int sign1 = 1;
    int sign2 = 1;
    QStringList halves = s.split("+");
    if(halves.length() == 1){
        halves = s.split("-", QString::SkipEmptyParts);
        sign2 = -1;
        if(s.at(0) == '-'){
            sign1 = -1;
        }
    }
    double r = halves.at(0).toDouble();
    QString tmp = halves.at(1);
    tmp.chop(1);
    double im = tmp.toDouble();
    return sign1*r + sign2*im*i;
}

void MainWindow::calc(){
    l = double(ui->lSlider->value())/2;
    ui->lLabel->setText(QString::number(l));
    int n = (l-1300)*2;
    if(n != 0){
        n--;
    }
    l = l*pow(10, -9);
    b = coeffs[n][0].real();
    ui->betaLabel->setText(QString::number(b));
    ui->b2Label->setText(QString::number(b/k));
    ga = sqrt(pow(b, 2) - pow(k*na, 2));
    kr = sqrt(pow(k*nr, 2) - pow(b, 2));
    gs = sqrt(pow(b, 2) - pow(k*ns, 2));
    gp = sqrt(pow(b, 2) - pow(k*np, 2));
    omega = k/sqrt(eps);
    if(chart->axes().contains(axisX)){
        chart->removeAxis(axisX);
    }
    if(chart->axes().contains(axisY)){
        chart->removeAxis(axisY);
    }
    axisY->setMin(-hr1 - 100*pow(10, -9));
    axisY->setMax(hr2 + hs + 100*pow(10, -9));
    if(chart->axes().contains(axisX)){
        chart->removeAxis(axisX);
    }
    if(chart->axes().contains(axisY)){
        chart->removeAxis(axisY);
    }
    chart->addAxis(axisX, Qt::AlignBottom);
    chart->addAxis(axisY, Qt::AlignLeft);
    QtCharts::QAbstractSeries *delTmp;
    foreach(delTmp, chart->series()){
        chart->removeSeries(delTmp);
    }
    int c = ui->pointSpinBox->value();
    qreal min = 100000000;
    qreal max = -10000000;
    double curr = -hr1 - 100*pow(10, -9);
    double step = (hr2 + hs + 100*pow(10, -9) -curr)/c;
    Hxi = new QtCharts::QLineSeries();
    Hxi->setName("Hx imag.");
    Hxr = new QtCharts::QLineSeries();
    Hxr->setName("Hx real");
    zero = new QtCharts::QLineSeries();
    zero->setName("zero");
    for(int it = 0; it < c; it++){
        zero->append(0, curr);
        if(curr <= -hr1){
            qreal tmp = (coeffs[n][2]*exp(gp*curr)).real();
            Hxr->append(tmp, curr);
            if(tmp > max){
                max = tmp;
            }else{
                if(tmp < min){
                    min = tmp;
                }
            }
            tmp = (coeffs[n][2]*exp(gp*curr)).imag();
            Hxi->append(tmp, curr);
            if(tmp > max){
                max = tmp;
            }else{
                if(tmp < min){
                    min = tmp;
                }
            }
        }else{
            if(curr <= 0){
                qreal tmp = (coeffs[n][3]*exp(i*curr*kr) + (coeffs[n][4])*exp(-i*curr*kr)).real();
                Hxr->append(tmp, curr);
                if(tmp > max){
                    max = tmp;
                }else{
                    if(tmp < min){
                        min = tmp;
                    }
                }
                tmp = (coeffs[n][3]*exp(i*curr*kr) + coeffs[n][4]*exp(-i*curr*kr)).imag();
                Hxi->append(tmp, curr);
                if(tmp > max){
                    max = tmp;
                }else{
                    if(tmp < min){
                        min = tmp;
                    }
                }
            }else{
                if(curr <= hs){
                    qreal tmp = (coeffs[n][5]*exp(curr*gs) + (coeffs[n][6])*exp(-curr*gs)).real();
                    Hxr->append(tmp, curr);
                    if(tmp > max){
                        max = tmp;
                    }else{
                        if(tmp < min){
                            min = tmp;
                        }
                    }
                    tmp = (coeffs[n][5]*exp(curr*gs) + (coeffs[n][6])*exp(-curr*gs)).imag();
                    Hxi->append(tmp, curr);
                    if(tmp > max){
                        max = tmp;
                    }else{
                        if(tmp < min){
                            min = tmp;
                        }
                    }
                }else{
                    if(curr <= hs + hr2){
                        qreal tmp = (coeffs[n][7]*exp(i*curr*kr) + (coeffs[n][8])*exp(-i*curr*kr)).real();
                        Hxr->append(tmp, curr);
                        if(tmp > max){
                            max = tmp;
                        }else{
                            if(tmp < min){
                                min = tmp;
                            }
                        }
                        tmp = (coeffs[n][7]*exp(i*curr*kr) + (coeffs[n][8])*exp(-i*curr*kr)).imag();
                        Hxi->append(tmp, curr);
                        if(tmp > max){
                            max = tmp;
                        }else{
                            if(tmp < min){
                                min = tmp;
                            }
                        }
                    }else{
                        qreal tmp = (coeffs[n][9]*exp(-ga*curr)).real();
                        Hxr->append(tmp, curr);
                        if(tmp > max && ui->hxBox->isChecked() && ui->riButton->isChecked()){
                            max = tmp;
                        }else{
                            if(tmp < min && ui->hxBox->isChecked() && ui->riButton->isChecked()){
                                min = tmp;
                            }
                        }
                        tmp = (coeffs[n][9]*exp(-ga*curr)).imag();
                        Hxi->append(tmp, curr);
                        if(tmp > max && ui->hxBox->isChecked() && ui->riButton->isChecked()){
                            max = tmp;
                        }else{
                            if(tmp < min && ui->hxBox->isChecked() && ui->riButton->isChecked()){
                                min = tmp;
                            }
                        }
                    }
                }
            }
        }
        curr += step;
    }
    zero->setColor(QColor(1, 1, 1));
    chart->addSeries(zero);
    zero->attachAxis(axisX);
    zero->attachAxis(axisY);
    if(ui->hxBox->isChecked() && ui->riButton->isChecked()){
        Hxi->setColor(QColor(255, 0, 0));
        chart->addSeries(Hxi);
        Hxi->attachAxis(axisX);
        Hxi->attachAxis(axisY);
        Hxr->setColor(QColor(150, 0, 0));
        chart->addSeries(Hxr);
        Hxr->attachAxis(axisX);
        Hxr->attachAxis(axisY);
    }else if (ui->hxBox->isChecked()) {
        Hxa = new QtCharts::QLineSeries();
        Hxa->setName("Hx arg.");
        Hxm = new QtCharts::QLineSeries();
        Hxm->setName("Hx abs");

        for(int it = 0; it< c; it++){
            double tmp;
            tmp = sqrt(pow(Hxi->at(it).x(), 2) + pow(Hxr->at(it).x(), 2));
            Hxm->append(tmp, Hxi->at(it).y());
            if(tmp > max){
                max = tmp;
            }else{
                if(tmp < min){
                    min = tmp;
                }
            }
            tmp = atan(Hxi->at(it).x()/Hxr->at(it).x());
            if(tmp > max){
                max = tmp;
            }else{
                if(tmp < min){
                    min = tmp;
                }
            }
            Hxa->append(tmp, Hxr->at(it).y());
        }

        Hxa->setColor(QColor(255, 0, 0));
        chart->addSeries(Hxa);
        Hxa->attachAxis(axisX);
        Hxa->attachAxis(axisY);
        Hxm->setColor(QColor(150, 0, 0));
        chart->addSeries(Hxm);
        Hxm->attachAxis(axisX);
        Hxm->attachAxis(axisY);
    }
    if(ui->eyBox->isChecked()){
        Eyi = new QtCharts::QLineSeries();
        Eyi->setName("Ey imag.");
        Eyr = new QtCharts::QLineSeries();
        Eyr->setName("Ey real");
        for(int it = 0; it< c; it++){
            double tmp;
            if(Hxi->at(it).y() <= -hr1){
                tmp = -b*Hxi->at(it).x()/(eps*omega*pow(np,2));
            }else{
                if((Hxi->at(it).y() <= 0) || ((Hxi->at(it).y() <= hr2 + hs) && ((Hxi->at(it).y() >= hs)))){
                    tmp = -b*Hxi->at(it).x()/(eps*omega*pow(nr,2));
                }else{
                    if(Hxi->at(it).y() <= hs){
                        tmp = -b*Hxi->at(it).x()/(eps*omega*pow(ns,2));
                    }else{
                        tmp = -b*Hxi->at(it).x()/(eps*omega*pow(na,2));
                    }
                }
            }
            if(tmp > max && ui->riButton->isChecked()){
                max = tmp;
            }else{
                if(tmp < min && ui->riButton->isChecked()){
                    min = tmp;
                }
            }
            Eyi->append(tmp, Hxi->at(it).y());

            if(Hxr->at(it).y() <= -hr1){
                tmp = -b*Hxr->at(it).x()/(eps*omega*pow(np,2));
            }else{
                if((Hxr->at(it).y() <= 0) || ((Hxr->at(it).y() <= hr2 + hs) && ((Hxr->at(it).y() >= hs)))){
                    tmp = -b*Hxr->at(it).x()/(eps*omega*pow(nr,2));
                }else{
                    if(Hxr->at(it).y() <= hs){
                        tmp = -b*Hxr->at(it).x()/(eps*omega*pow(ns,2));
                    }else{
                        tmp = -b*Hxr->at(it).x()/(eps*omega*pow(na,2));
                    }
                }
            }
            if(tmp > max && ui->riButton->isChecked()){
                max = tmp;
            }else{
                if(tmp < min && ui->riButton->isChecked()){
                    min = tmp;
                }
            }
            Eyr->append(tmp, Hxr->at(it).y());
        }
        if(ui->riButton->isChecked()){
            Eyi->setColor(QColor(0, 255, 0));
            chart->addSeries(Eyi);
            Eyi->attachAxis(axisX);
            Eyi->attachAxis(axisY);
            Eyr->setColor(QColor(0, 150, 0));
            chart->addSeries(Eyr);
            Eyr->attachAxis(axisX);
            Eyr->attachAxis(axisY);
        }else if (ui->eyBox->isChecked()) {
            Eya = new QtCharts::QLineSeries();
            Eya->setName("Ey arg.");
            Eym = new QtCharts::QLineSeries();
            Eym->setName("Ey abs");

            for(int it = 0; it< c; it++){
                double tmp;
                tmp = sqrt(pow(Eyi->at(it).x(), 2) + pow(Eyr->at(it).x(), 2));
                Eym->append(tmp, Eyi->at(it).y());
                if(tmp > max){
                    max = tmp;
                }else{
                    if(tmp < min){
                        min = tmp;
                    }
                }
                tmp = atan(Eyi->at(it).x()/Eyr->at(it).x());
                if(tmp > max){
                    max = tmp;
                }else{
                    if(tmp < min){
                        min = tmp;
                    }
                }
                Eya->append(tmp, Eyr->at(it).y());
            }

            Eya->setColor(QColor(0, 255, 0));
            chart->addSeries(Eya);
            Eya->attachAxis(axisX);
            Eya->attachAxis(axisY);
            Eym->setColor(QColor(0, 150, 0));
            chart->addSeries(Eym);
            Eym->attachAxis(axisX);
            Eym->attachAxis(axisY);
        }
    }
    if(ui->ezBox->isChecked()){
        Ezi = new QtCharts::QLineSeries();
        Ezi->setName("Ez imag.");
        Ezr = new QtCharts::QLineSeries();
        Ezr->setName("Ez real");

        for(int it = 0; it< c; it++){
            std::complex<double> tmp;
            if(Hxi->at(it).y() <= -hr1){
                tmp = (gp*coeffs[n][2]*exp(gp*Hxi->at(it).y()))/pow(np,2);
            }else{
                if(Hxi->at(it).y() <= 0){
                    tmp = (i*kr*(coeffs[n][3]*exp(i*Hxi->at(it).y()*kr) - coeffs[n][4]*exp(-i*Hxi->at(it).y()*kr)))/pow(nr,2);
                }else{
                    if(Hxi->at(it).y() <= hs){
                        tmp = (gs*(coeffs[n][5]*exp(Hxi->at(it).y()*gs) - (coeffs[n][6])*exp(-Hxi->at(it).y()*gs)))/pow(ns,2);
                    }else{
                        if(Hxi->at(it).y() <= hs+hr2){
                            tmp = (i*kr*(coeffs[n][7]*exp(i*Hxi->at(it).y()*kr) - coeffs[n][8]*exp(-i*Hxi->at(it).y()*kr)))/pow(nr,2);
                        }else{
                            tmp = -ga*(coeffs[n][9]*exp(-ga*Hxi->at(it).y()))/pow(na,2);
                        }
                    }
                }
            }
            tmp = tmp/(i*eps*omega);
            if(tmp.real() > max && ui->riButton->isChecked()){
                max = tmp.real();
            }else{
                if(tmp.real() < min && ui->riButton->isChecked()){
                    min = tmp.real();
                }
            }
            Ezr->append(tmp.real(), Hxr->at(it).y());
            if(tmp.imag() > max && ui->riButton->isChecked()){
                max = tmp.imag();
            }else{
                if(tmp.imag() < min && ui->riButton->isChecked()){
                    min = tmp.imag();
                }
            }
            Ezi->append(tmp.imag(), Hxi->at(it).y());
        }

        if(ui->riButton->isChecked()){
            Ezi->setColor(QColor(0, 0, 255));
            chart->addSeries(Ezi);
            Ezi->attachAxis(axisX);
            Ezi->attachAxis(axisY);
            Ezr->setColor(QColor(0, 0, 150));
            chart->addSeries(Ezr);
            Ezr->attachAxis(axisX);
            Ezr->attachAxis(axisY);
        }else if (ui->ezBox->isChecked()) {
            Eza = new QtCharts::QLineSeries();
            Eza->setName("Ey arg.");
            Ezm = new QtCharts::QLineSeries();
            Ezm->setName("Ey abs");

            for(int it = 0; it< c; it++){
                double tmp;
                tmp = sqrt(pow(Ezi->at(it).x(), 2) + pow(Ezr->at(it).x(), 2));
                Ezm->append(tmp, Ezi->at(it).y());
                if(tmp > max){
                    max = tmp;
                }else{
                    if(tmp < min){
                        min = tmp;
                    }
                }
                tmp = atan(Ezi->at(it).x()/Ezr->at(it).x());
                if(tmp > max){
                    max = tmp;
                }else{
                    if(tmp < min){
                        min = tmp;
                    }
                }
                Eza->append(tmp, Ezr->at(it).y());
            }

            Eza->setColor(QColor(0, 0, 255));
            chart->addSeries(Eza);
            Eza->attachAxis(axisX);
            Eza->attachAxis(axisY);
            Ezm->setColor(QColor(0, 0, 150));
            chart->addSeries(Ezm);
            Ezm->attachAxis(axisX);
            Ezm->attachAxis(axisY);
        }
    }

    axisX->setMin(min);
    axisX->setMax(max);
}

void MainWindow::save(){
    outFile = new QFile(ui->lineEdit->text());
    outFile->open(QIODevice::WriteOnly);
    outStream = new QTextStream(outFile);
    int c = ui->pointSpinBox->value();
    double curr = -hr1 - 100*pow(10, -9);
    double step = (hr2 + hs + 100*pow(10, -9) -curr)/c;
    *outStream << "wawelength = " << ui->lLabel->text() << endl;
    *outStream << qSetFieldWidth(14) << left;
    *outStream << "y";
    if(ui->hxBox->isChecked()){
        if(ui->riButton->isChecked()){
            *outStream << "Re(Hx)";
            *outStream << "Im(Hx)";
        }else{
            *outStream << "|Hx|";
            *outStream << "arg(Hx)";
        }
    }
    if(ui->eyBox->isChecked()){
        if(ui->riButton->isChecked()){
            *outStream << "Re(Ey)";
            *outStream << "Im(Ey)";
        }else{
            *outStream << "|Ey|";
            *outStream << "arg(Ey)";
        }
    }
    if(ui->ezBox->isChecked()){
        if(ui->riButton->isChecked()){
            *outStream << "Re(Ez)";
            *outStream << "Im(Ez)";
        }else{
            *outStream << "|Ez|";
            *outStream << "arg(Ez)";
        }
    }
    for(int it = 0; it < c; it++){
        *outStream << qSetFieldWidth(1) << "\n";
        *outStream << qSetFieldWidth(14) << left;
        *outStream << curr;
        if(ui->hxBox->isChecked()){
            if(ui->riButton->isChecked()){
                *outStream << Hxr->at(it).x();
                *outStream << Hxi->at(it).x();
            }else{
                *outStream << Hxm->at(it).x();
                *outStream << Hxa->at(it).x();
            }
        }
        if(ui->eyBox->isChecked()){
            if(ui->riButton->isChecked()){
                *outStream << Eyr->at(it).x();
                *outStream << Eyi->at(it).x();
            }else{
                *outStream << Eym->at(it).x();
                *outStream << Eya->at(it).x();
            }
        }
        if(ui->ezBox->isChecked()){
            if(ui->riButton->isChecked()){
                *outStream << Ezr->at(it).x();
                *outStream << Ezi->at(it).x();
            }else{
                *outStream << Ezm->at(it).x();
                *outStream << Eza->at(it).x();
            }
        }
        curr += step;
    }
    outFile->close();
}

void MainWindow::deb(std::complex<double> n){
    qDebug() << n.real() << "+ i*" << n.imag();
}

void MainWindow::resizeEvent(QResizeEvent* event){
    QMainWindow::resizeEvent(event);
    chartView->resize(ui->chartWidget->size());
}

<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>
<script type='text/x-mathjax-config'>
	MathJax.Hub.Config({
		tex2jax: {
			inlineMath: [['$', '$']],
			dsplayMath: [ ['$$','$$'],["\\[","\\]"] ]
		}
	});
</script>

# ActFEMtide マニュアル
## 潮汐起源電磁場のための有限要素法シミュレーションコード

#### 2026年7月23日現在

南 拓人
tminami@port.kobe-u.ac.jp
（神戸大学大学院理学研究科）

### 目次

1. [はじめに](#Introduction)
2. [重要事項](#ImportantNotes)
3. [必要な環境](#Required_environments)
4. [サンプルシミュレーションの実行](#Run_sample_simulations)  
 - サンプル1: [鬼海カルデラ (Kikai)](#Tohoku_small)
5. [支配方程式と離散化](#Governing_equation_and_discretization)
   - [5.1 支配方程式](#sec5_1)
   - [5.2 弱形式](#sec5_2)
   - [5.3 辺要素による離散化](#sec5_3)
   - [5.4 潮汐起電力項（右辺ベクトル）](#sec5_4)
   - [5.5 スケーリングと単位](#sec5_5)
   - [5.6 電場・磁場・電流の計算](#sec5_6)

## 1 はじめに <a id="Introduction"></a>
__ActFEMtide__ は、海洋中の起電力を介して潮汐起源の電磁場を計算するためのシミュレーションコードである。メインのソルバーは __ActFEM__ [(Minami et al. 2018)](#Minami2018) をベースとしており、励起源を電流ダイポールから潮汐起源の起電力（$\mathbf{v}\times \mathbf{F}$）に置き換えたものである。
__ActFEMtide__ は Fortran で書かれており、疎行列直接解法ソルバー PARDISO を使用するために Intel Math Kernel Library (mkl) を利用する。ActFEMtide は PARDISO による並列計算に openMP を、複数周波数に対する並列計算に MPI を使用できる。
__ActFEMtide__ は __TMTGEM__ [(Minami et al. 2017)](#Minami2017) のために開発されたメッシュ生成部を使用している。
ActFEMtide は 1: メッシュ生成部、2: シミュレーション部の2つの部分から構成される。__ActFEMtide__ は現在、[(Egbert and Erofeeva, 2002)](#EE2002) による潮汐モデルの出力を使用することを前提としている。
__ActFEMtide__ はメッシュ生成に [Gmsh](https://www.soest.hawaii.edu/gmt/) を使用し、[Oishi et al. (2013)](#O2013) に記載された extrude アルゴリズムを用いる。これはもともと __TMTGEM__ のために開発されたものである。extrusion のためのサブルーチンの多くは、流体シミュレーションのオープンソースコードである [Fluidity](http://fluidityproject.github.io/) に由来する。

なお、本マニュアルには __TMTGEM__ マニュアル（v1.3）からの流用箇所が多く含まれることに注意されたい。

## 2 重要事項 <a id="ImportantNotes"></a>
座標系: X: 東向き, Y: 北向き, Z: 上向き

__支配方程式:__ [5 支配方程式と離散化](#Governing_equation_and_discretization) を参照（[Minami et al. (2018, EPS)](#Minami2018) の送信源項を起電力項 $\sigma\,\mathbf{v}\times \mathbf{F}$ に変更）

__入力する地形データ:__ 経度 [deg], 緯度 [deg], 高度 [m, __下向き正__]

__使用する潮汐モデル:__ [Egbert and Erofeeva (2002)](#EE2002)

## 3 必要な環境 <a id="Required_environments"></a>
以下のパッケージをPCにインストールしておくこと:

- __Gmsh__ [(http://gmsh.info/)](http://gmsh.info/) （メッシュ生成用）

- __mkl ライブラリ付きの Intel fortran コンパイラ__ （“ifort –mkl=parallel ***” が使用できることを確認すること）

- __GMT__ ([generic mapping tool](https://www.soest.hawaii.edu/gmt/)) （シミュレーション結果の可視化用）
ghostscript （TMTGEM/Tohoku/flow 内の plot_z.sh などで使用する “gv” コマンドのみに必要）

Debian 系の Linux ディストリビューションが使用できる場合、"sudo apt install ***" によって上記の環境をすべてインストールできる。

## 4 "Kikai/" 内のサンプルコードの実行 <a id="Run_sample_simulations"></a>
とにかく ActFEMtide を実行してみよう！！
ActFEMtide のホームフォルダ内に、サンプルシミュレーション用の Kikai/ フォルダがある。
そこにあるサンプルコードは、以下のいくつかのステップで実行される。

### Step0 潮汐モデルの準備
[(Egbert and Erofeeva, 2002)](#EE2002) のコードを ActFEMtide/mkfvxyz/ 内でコンパイルする
####  
    $cd mkfvxyz
    $cd OTPS
    $make           (潮汐モデルを使用するためのコードをコンパイルする)

triton で作業している場合は、ActFEMtide/mkfvxyz/OTPS/ 内で以下のコマンドによりデータをコピーする
####
    $./wget.sh
DATA フォルダ内に、以下のファイルが準備されていることを確認すること:


### Step 1: 「地形データの準備」(ActFEMtide/Kikai/topo/)
#### etopo の grd ファイルを ascii 形式の *.xyz ファイルに変換する
    $cd Kikai/topo
    $./mk_etopo_kikai.sh        (etopo_kikai-l.xyz が生成される)

### Step 2: 四面体メッシュの生成 (ActFEMtide/Kikai/mesh/)
#### 四面体メッシュの生成
    $cd ../mesh
    $./tetmeshgen.sh             (em3d.msh 等が生成される)
使用しているPCのスペックによっては、これにはしばらく時間がかかる。多くのファイルが生成されるが、鍵となる三次元四面体メッシュファイルは "em3d.msh" である。メッシュファイルをローカルにダウンロードし、以下のステップで確認すること。

### Step 3: 生成されたメッシュファイルの確認
#### ActFEM/Kikai/mesh 内のファイルをダウンロードした後
    $gmsh em3d.msh

メッシュの生成に成功していれば、メッシュを見ることができる。

## Step 4: 潮汐流と背景磁場の準備
####
    $cd ../fvxyz
    $./mkfvxyz_mesh.sh

## Step 5: ActFEMtide の実行
####
    $cd ../fwd
    $./run_fwd.sh

## Step 6: 結果の確認
    $./plot_bxyz.sh
    $./plot_ixyh.sh
以下の2つの図が生成される。

![Fig.1](./images/bxyz.png)
図1. M2潮汐による潮汐起源磁場。3成分すべて海底における値である。この図は plot_bxyz.sh により生成される。

![Fig.2](./images/ixyh.png)
図2. M2潮汐による磁場・電場の実部および虚部。この図は plot_ixyh.sh により生成される。

## 5 支配方程式と離散化 <a id="Governing_equation_and_discretization"></a>
ActFEMtide のソルバーは __ActFEM__ [(Minami et al. 2018)](#Minami2018) と同じく、ベクトルポテンシャル $\mathbf{A}$ に関する誘導方程式を、四面体メッシュ上の辺要素（Nédélec 一次要素）有限要素法で周波数領域において解く。ActFEM との違いは右辺の励起源のみであり、ActFEM の送信源電流 $\mathbf{i}_s$ の代わりに、海水の運動により生じる起電力 $\sigma\,\mathbf{v}\times\mathbf{F}$ を与える。以下、離散化の導出は ActFEM マニュアル（v1.5、付録A・B）に従う。

### 5.1 支配方程式 <a id="sec5_1"></a>
導電率 $\sigma$ の海水が、背景磁場（地球主磁場）$\mathbf{F}$ の中を速度 $\mathbf{v}$ で運動するとき、オームの法則は

$$\mathbf{j} = \sigma\left( \mathbf{E} + \mathbf{v}\times\mathbf{F} \right)$$

となる。変位電流を無視したアンペールの法則 $\nabla\times\mathbf{B} = \mu\mathbf{j}$ に $\mathbf{B} = \nabla\times\mathbf{A}$ を代入し、ActFEM と同じ規約 $\mathbf{E} = i\omega\mathbf{A}$ で電場を表すと、ActFEMtide が解く支配方程式

$$\nabla\times\left( \nabla\times\mathbf{A} \right) + i\omega\mu\sigma\mathbf{A} = \mu\sigma\left( \mathbf{v}\times\mathbf{F} \right) \cdot \cdot \cdot (5.1)$$

が得られる。ここで $\omega$ は潮汐分潮の角周波数、$\mu$ は透磁率（真空の透磁率 $\mu_0$ を使用）である。ActFEM の式 $\nabla\times\nabla\times\mathbf{A} + i\omega\mu\sigma\mathbf{A} = \mu\mathbf{i}_s$ と比べると、送信源電流 $\mathbf{i}_s$ が起電力による電流 $\sigma\,\mathbf{v}\times\mathbf{F}$ に置き換わっている。

- $\mathbf{v}$ は潮汐流速であり、[Egbert and Erofeeva (2002)](#EE2002) の潮汐モデルから得た複素振幅（振幅と位相）で与える。水平2成分のみを使用し、鉛直成分は $v_z = 0$ とする。
- $\mathbf{F}$ は背景磁場（3成分、実数）である。
- 起電力項は海洋要素（物理グループ 2）にのみ存在する。

### 5.2 弱形式 <a id="sec5_2"></a>
ガラーキン法では、試験関数 $\mathbf{w}$（辺要素基底関数）をかけて計算領域 $\Omega$ で積分した

$$\int_{\Omega}\mathbf{w}\cdot\left\{ \nabla\times\nabla\times\mathbf{A} + i\omega\mu\sigma\mathbf{A} - \mu\sigma\left( \mathbf{v}\times\mathbf{F} \right) \right\}dV = 0 \cdot \cdot \cdot (5.2)$$

を解く。第1項を部分積分すると、

$$\int_{\Omega}\mathbf{w}\cdot\left( \nabla\times\nabla\times\mathbf{A} \right)dV = \int_{\partial\Omega}\left( \mathbf{B}\times\mathbf{w} \right)\cdot\mathbf{n}\,dS + \int_{\Omega}\left( \nabla\times\mathbf{w} \right)\cdot\left( \nabla\times\mathbf{A} \right)dV$$

となる。計算領域の外側境界 $\partial\Omega$ では ActFEM と同様に $\mathbf{n}\cdot\left( \mathbf{B}\times\mathbf{w} \right) = \mathbf{w}\cdot\left( \mathbf{n}\times\mathbf{B} \right) = 0$ とおけるため（境界を貫く電流がないことに対応）、境界積分は消える。したがって弱形式は

$$\int_{\Omega}\left( \nabla\times\mathbf{w} \right)\cdot\left( \nabla\times\mathbf{A} \right)dV + i\omega\mu\int_{\Omega}\sigma\,\mathbf{w}\cdot\mathbf{A}\,dV = \mu\int_{\Omega}\sigma\,\mathbf{w}\cdot\left( \mathbf{v}\times\mathbf{F} \right)dV \cdot \cdot \cdot (5.3)$$

となる。境界では $\mathbf{A}$ の接線成分を 0 とする（Dirichlet 条件）。計算領域は十分大きくとる必要がある。

### 5.3 辺要素による離散化 <a id="sec5_3"></a>
計算領域を四面体要素に分割し、$e$ 番目の要素 $\Omega_e$ 内のベクトルポテンシャルを、6本の辺に対応する辺要素基底関数 $\mathbf{w}_j$ で

$$\mathbf{A}_e = \sum_{j=1}^{6}\mathbf{w}_j\left\lbrack Al \right\rbrack_j, \qquad \mathbf{w}_j = \mathbf{w}_{kl} = \lambda_k\nabla\lambda_l - \lambda_l\nabla\lambda_k \cdot \cdot \cdot (5.4)$$

と展開する。ここで $\lambda_k$ は節点 $k$ の一次節点基底関数（体積座標）、未知数 $\left\lbrack Al \right\rbrack_j$ は辺 $j$ に沿った $\mathbf{A}$ の線積分である。辺番号 $j$ と節点番号 $(k, l)$、対辺 $(m, n)$ の対応は下図と表5.1の通りであり、ActFEM と同じである。

<img src="./images/edge_numbering.png" width="55%">

表5.1.

| 辺番号 $j$ | $k$ | $l$ | $m$ | $n$ |
|:----------:|:---:|:---:|:---:|:---:|
| 1          | 1   | 2   | 3   | 4   |
| 2          | 2   | 3   | 1   | 4   |
| 3          | 1   | 3   | 4   | 2   |
| 4          | 1   | 4   | 2   | 3   |
| 5          | 2   | 4   | 3   | 1   |
| 6          | 3   | 4   | 1   | 2   |

式(5.4)を式(5.3)に代入すると、要素ごとの方程式

$$\left\{ \left\lbrack M_e \right\rbrack + i\omega\mu\sigma_e\left\lbrack N_e \right\rbrack \right\}\left\lbrack \mathbf{Al}_e \right\rbrack = \mathbf{s}_e \cdot \cdot \cdot (5.5)$$

が得られる。$\sigma_e$ は要素 $e$ の導電率（要素内一定）である。要素行列は解析的に計算できる。

$$\left\lbrack M_e \right\rbrack_{ij} = \int_{\Omega_e}\left( \nabla\times\mathbf{w}_i \right)\cdot\left( \nabla\times\mathbf{w}_j \right)dV = \frac{1}{9v}\,\mathbf{x}_{mn}\cdot\mathbf{x}_{m'n'} \cdot \cdot \cdot (5.6)$$

$$\left\lbrack N_e \right\rbrack_{ij} = \int_{\Omega_e}\mathbf{w}_i\cdot\mathbf{w}_j\,dV = \nabla\lambda_l\cdot\nabla\lambda_{l'}\int_{\Omega_e}\lambda_k\lambda_{k'}dV - \nabla\lambda_l\cdot\nabla\lambda_{k'}\int_{\Omega_e}\lambda_k\lambda_{l'}dV - \nabla\lambda_k\cdot\nabla\lambda_{l'}\int_{\Omega_e}\lambda_l\lambda_{k'}dV + \nabla\lambda_k\cdot\nabla\lambda_{k'}\int_{\Omega_e}\lambda_l\lambda_{l'}dV \cdot \cdot \cdot (5.7)$$

ここで $(k,l,m,n)$ は辺 $i$、$(k',l',m',n')$ は辺 $j$ に対応する節点番号、$v$ は四面体の体積、$\mathbf{x}_{mn}$ は節点 $m$ から $n$ へのベクトルである。導出には次の関係を用いる。

$$\nabla\times\mathbf{w}_j = 2\nabla\lambda_k\times\nabla\lambda_l = \frac{1}{3v}\mathbf{x}_{mn}, \qquad \int_{\Omega_e}\lambda_k\lambda_l\,dV = \begin{cases} v/20 & (k\neq l) \\ v/10 & (k = l) \end{cases} \cdot \cdot \cdot (5.8)$$

辺の向きは、要素ローカルの辺と全体の辺（line）の向きが逆の場合に符号 $-1$ をかけて揃える。要素行列を全要素について重ね合わせて得られる全体の連立一次方程式（複素対称）を、PARDISO で解く。

各要素の導電率 $\sigma_e$ は物理グループで決まる。

| 物理グループ | 領域 | $\sigma_e$ |
|:----:|:----|:----|
| 1 | 空気 | 制御ファイルで与える $\sigma_\mathrm{air}$ |
| 2 | 海洋 | 3.33 S/m（ソース内で固定） |
| 3 以上 | 陸・海底下 | 制御ファイルの値、または比抵抗ファイル（cond.msh）の値 |

### 5.4 潮汐起電力項（右辺ベクトル） <a id="sec5_4"></a>
右辺ベクトルの要素寄与は

$$s_i^e = \mu\sigma_e\int_{\Omega_e}\mathbf{w}_i\cdot\left( \mathbf{v}\times\mathbf{F} \right)dV \cdot \cdot \cdot (5.9)$$

である。ActFEM では送信源が線電流（ヘヴィサイド関数）であったため、この積分は送信ケーブルが貫く要素だけで評価していた。ActFEMtide では起電力が海洋全体に分布するため、海洋要素すべてで評価する。各節点 $p$ における流速 $\mathbf{v}_p$ と背景磁場 $\mathbf{F}_p$ から節点値 $\left( \mathbf{v}\times\mathbf{F} \right)_p = \mathbf{v}_p\times\mathbf{F}_p$ を計算し、要素内では節点基底関数で一次補間する。

$$\left( \mathbf{v}\times\mathbf{F} \right)(\mathbf{r}) \simeq \sum_{p=1}^{4}\lambda_p(\mathbf{r})\left( \mathbf{v}\times\mathbf{F} \right)_p \cdot \cdot \cdot (5.10)$$

式(5.10)と $\mathbf{w}_i = \lambda_k\nabla\lambda_l - \lambda_l\nabla\lambda_k$ を式(5.9)に代入すると、$\nabla\lambda$ が要素内で一定であることから

$$s_i^e = \mu\sigma_e\sum_{p=1}^{4}\left\{ \left( \int_{\Omega_e}\lambda_k\lambda_p\,dV \right)\nabla\lambda_l\cdot\left( \mathbf{v}\times\mathbf{F} \right)_p - \left( \int_{\Omega_e}\lambda_l\lambda_p\,dV \right)\nabla\lambda_k\cdot\left( \mathbf{v}\times\mathbf{F} \right)_p \right\} \cdot \cdot \cdot (5.11)$$

となり、式(5.8)の積分公式を用いて解析的に計算できる。$\mathbf{v}$ が複素数のため、$s_i^e$ も複素数となる。海洋以外の要素では $s_i^e = 0$ である。

### 5.5 スケーリングと単位 <a id="sec5_5"></a>
ActFEM と同様に、スケーリング長 $L = 1000$ m を用いて座標を km 単位で扱う（ActFEM マニュアル付録B）。$\nabla = \nabla'/L$、$\mathbf{w} = \mathbf{w}'/L$、$dV = L^3dV'$ と置き換えると、式(5.5)は

$$\sum_{j=1}^{6}\left\{ \left\lbrack M_e' \right\rbrack_{ij} + L^2 i\omega\mu\sigma_e\left\lbrack N_e' \right\rbrack_{ij} \right\}\left\lbrack Al' \right\rbrack_j = \mu\sigma_e L^2\int_{\Omega_e}\mathbf{w}_i'\cdot\left( \mathbf{v}\times\mathbf{F} \right)dV' \cdot \cdot \cdot (5.12)$$

となる。ここで $'$ は km 単位の座標で構成した量を表す。

入力の単位は流速 $\mathbf{v}$ [m/s]、背景磁場 $\mathbf{F}$ [nT] であり、$\mathbf{v}\times\mathbf{F}$ は [nV/m] で与えられる。一方、ActFEMtide は電場を [mV/km]（= [μV/m]）で出力する。このためコード（src/solver_mpi/forward_bxyz_mpi.f90 の GENMAT）では、式(5.12)の右辺に次の係数をかけている。

$$L^2 \times 10^{6} \times 10^{-9} = 10^{6}\,[\mathrm{m^2/km^2}] \times 10^{6}\,[(\mathrm{mV/km})/(\mathrm{V/m})] \times 10^{-9}\,[(\mathrm{V/m})/(\mathrm{nV/m})] = 10^{3}$$

### 5.6 電場・磁場・電流の計算 <a id="sec5_6"></a>
連立方程式を解いて $\left\lbrack Al \right\rbrack$ が求まった後、観測点を含む要素内で式(5.4)により $\mathbf{A}$ を補間し、電場を $\mathbf{E} = i\omega\mathbf{A}$、磁場を $\mathbf{B} = \nabla\times\mathbf{A}$ から求める。

海洋中の電流には、誘導電場 $\mathbf{E}$ による寄与と起電力 $\mathbf{v}\times\mathbf{F}$ による寄与の両方が含まれる（5.1節のオームの法則）。plot_ixyh.sh で描く水平電流（深さ積分）は、水深 $h$ を用いて

$$\mathbf{I}_h = \sigma h\left( \mathbf{E}_h + \left( \mathbf{v}\times\mathbf{F} \right)_h \right)$$

として計算し、$\sigma h\mathbf{E}_h$ と $\sigma h(\mathbf{v}\times\mathbf{F})_h$ の各成分も出力する（src/solver_mpi/n_ebfem_bxyz_mpi.f90）。

## 参考文献
- Egbert, G. D., & Erofeeva, S. Y. (2002). Efficient inverse modeling of barotropic ocean tides. Journal of Atmospheric and Oceanic technology, 19(2), 183-204.<a id="EE2022"></a>
- Minami, T., Utsugi, M., Utada, H., Kagiyama, T., & Inoue, H. (2018). Temporal variation in the resistivity structure of the first Nakadake crater, Aso volcano, Japan, during the magmatic eruptions from November 2014 to May 2015, as inferred by the ACTIVE electromagnetic monitoring system. Earth, Planets and Space, 70(1), 138. <a id="Minami2018"></a>
- Minami, T., Toh, H., Ichihara, H., & Kawashima, I. (2017). Three‐Dimensional Time Domain Simulation of Tsunami‐Generated Electromagnetic Fields: Application to the 2011 Tohoku Earthquake Tsunami. Journal of Geophysical Research: Solid Earth, 122(12), 9559-9579. <a id="M2017"></a>

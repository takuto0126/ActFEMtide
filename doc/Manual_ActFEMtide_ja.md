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
 - サンプル1: [喜界島 (Kikai)](#Tohoku_small)

## 1 はじめに <a id="Introduction"></a>
__ActFEMtide__ は、海洋中の起電力を介して潮汐起源の電磁場を計算するためのシミュレーションコードである。メインのソルバーは __ActFEM__ [(Minami et al. 2018)](#Minami2018) をベースとしており、電気ダイポール源の項を潮汐起源の起電力（$\mathbf{v}\times \mathbf{F}$）に置き換えたものである。
__ActFEMtide__ は Fortran で書かれており、疎行列直接解法ソルバー PARDISO を使用するために Intel Math Kernel Library (mkl) を利用する。ActFEMtide は PARDISO による並列計算に openMP を、複数周波数に対する並列計算に MPI を使用できる。
__ActFEMtide__ は __TMTGEM__ [(Minami et al. 2017)](#Minami2017) のために開発されたメッシュ生成部を使用している。
ActFEMtide は 1: メッシュ生成部、2: シミュレーション部の2つの部分から構成される。__ActFEMtide__ は現在、[(Egbert and Erofeeva, 2002)](#EE2002) による潮汐モデルの出力を使用することを前提としている。サンプルコードには、[Kawashima and Toh (2016)](#KT2016) により改良された COMCOT ver 1.7 が含まれている。
__ActFEMtide__ はメッシュ生成に [Gmsh](https://www.soest.hawaii.edu/gmt/) を使用し、[Oishi et al. (2013)](#O2013) に記載された extrude アルゴリズムを用いる。これはもともと __TMTGEM__ のために開発されたものである。extrusion のためのサブルーチンの多くは、流体シミュレーションのオープンソースコードである [Fluidity](http://fluidityproject.github.io/) に由来する。

なお、本マニュアルには __TMTGEM__ マニュアル（v1.3）からの流用箇所が多く含まれることに注意されたい。

## 2 重要事項 <a id="ImportantNotes"></a>
座標系: X: 東向き, Y: 北向き, Z: 上向き

__支配方程式:__ [Minami et al. (2018, EPS)](#Minami2018) を参照（起電力項を $\mathbf{v}\times \mathbf{F}$ に変更）

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

## 参考文献
- Egbert, G. D., & Erofeeva, S. Y. (2002). Efficient inverse modeling of barotropic ocean tides. Journal of Atmospheric and Oceanic technology, 19(2), 183-204.<a id="EE2022"></a>
- Minami, T., Utsugi, M., Utada, H., Kagiyama, T., & Inoue, H. (2018). Temporal variation in the resistivity structure of the first Nakadake crater, Aso volcano, Japan, during the magmatic eruptions from November 2014 to May 2015, as inferred by the ACTIVE electromagnetic monitoring system. Earth, Planets and Space, 70(1), 138. <a id="Minami2018"></a>
- Minami, T., Toh, H., Ichihara, H., & Kawashima, I. (2017). Three‐Dimensional Time Domain Simulation of Tsunami‐Generated Electromagnetic Fields: Application to the 2011 Tohoku Earthquake Tsunami. Journal of Geophysical Research: Solid Earth, 122(12), 9559-9579. <a id="M2017"></a>

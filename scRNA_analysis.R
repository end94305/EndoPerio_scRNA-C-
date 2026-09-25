library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
# ===================================================================
# 1. データの準備: 根尖性歯周炎 (Apical) のデータのみを抽出
# ===================================================================
# ※ disease_status 列で Apical_Periodontitis を抽出します
seurat_apical <- subset(seurat_combined, disease_status == "Apical_Periodontitis")
# ※重要: CellChatでは細胞のアイデンティティ（Macrophage, Fibroblast等）が
# Identsにセットされている必要があります。もし別列（例: cell_type）にある場合は
# 以下のアクティブ化を行ってください（シャープを外して実行）。
# Idents(seurat_apical) <- seurat_apical$cell_type
# ===================================================================
# 2. CellChatパイプラインの実行（計算に少し時間がかかります）
# ===================================================================
# オブジェクトの作成
cellchat_apical <- createCellChat(object = seurat_apical, group.by = "ident")
# 1. 念のため、現在アクティブな細胞型(Idents)をメタデータに新しい列として明示的に追加します
seurat_apical$cell_type <- Idents(seurat_apical)
# 2. Seurat V5の正しい文法(layer)を使って、正規化済みデータとメタデータを個別に抽出します
data.input <- GetAssayData(seurat_apical, assay = "RNA", layer = "data")
meta_data <- seurat_apical@meta.data
# 3. 抽出したマトリックスとメタデータを使ってCellChatオブジェクトを作成します
cellchat_apical <- createCellChat(object = data.input, meta = meta_data, group.by = "cell_type")
# 4. データベースのセットアップ（ヒト）
CellChatDB <- CellChatDB.human
cellchat_apical@DB <- CellChatDB
# 5. 前処理と過剰発現遺伝子の特定
cellchat_apical <- subsetData(cellchat_apical)
cellchat_apical <- identifyOverExpressedGenes(cellchat_apical)
cellchat_apical <- identifyOverExpressedInteractions(cellchat_apical)
# 6. 通信確率の計算 (triMean法)
cellchat_apical <- computeCommunProb(cellchat_apical, type = "triMean")
colnames(seurat_apical@meta.data)
# 1. 正しい細胞ラベルの列（cellchat_labels）を指定してCellChatオブジェクトを作成
cellchat_apical <- createCellChat(object = data.input, meta = meta_data, group.by = "cellchat_labels")
# ※ここでコンソールに「The cell groups used for CellChat analysis are ...」と表示され、
# Macrophage, Fibroblast, Epithelial 等の細胞名が複数並んでいることを確認してください。
# 2. データベースのセットアップからパイプラインを再開
CellChatDB <- CellChatDB.human
cellchat_apical@DB <- CellChatDB
cellchat_apical <- subsetData(cellchat_apical)
cellchat_apical <- identifyOverExpressedGenes(cellchat_apical)
# 必要なパッケージの読み込み
library(CellChat)
library(ggplot2)
# 論文で言及した特異的なシグナル経路を指定
target_pathways <- c("COLLAGEN", "FN1", "LAMININ", "PTN", "APP", "MIF")
# すでに環境内にある cellchat_apical を使用してバブルプロットを作成
p_bubble <- netVisual_bubble(
cellchat_apical,
sources.use = "Macrophage",
targets.use = c("Fibroblast", "Epithelial"),
signaling = target_pathways,
remove.isolate = FALSE,
title.name = "Signaling from Macrophages in Apical Periodontitis"
) +
theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
# 必要なパッケージの読み込み
library(CellChat)
library(ggplot2)
# 論文で言及した特異的なシグナル経路を指定
target_pathways <- c("COLLAGEN", "FN1", "LAMININ", "PTN", "APP", "MIF")
# すでに環境内にある cellchat_apical を使用してバブルプロットを作成
p_bubble <- netVisual_bubble(
cellchat_apical,
sources.use = "Macrophage",
targets.use = c("Fibroblast", "Epithelial"),
signaling = target_pathways,
remove.isolate = FALSE,
title.name = "Signaling from Macrophages in Apical Periodontitis"
) +
theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
# 1. cellchat_apical に入っている細胞名の確認
levels(cellchat_apical@idents$joint)
# 1. 必要なパッケージの読み込み
library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)
# 2. 根尖性歯周炎 (Apical_Periodontitis) のデータだけを完全に独立したSeuratオブジェクトとして抽出
seurat_apical <- subset(seurat_combined, disease_status == "Apical_Periodontitis")
# 3. CellChat用の細胞ラベルを明確にするため、Identsに適切な列（例: cellchat_labels や seurat_clusters）をセット
# ※もしエラーが出る場合は、meta.dataにある細胞種が入っている列名に変更してください
Idents(seurat_apical) <- seurat_apical$cellchat_labels
# 4. CellChatオブジェクトの作成（マトリックスとメタデータを直接渡す安全な方法）
data.input <- GetAssayData(seurat_apical, assay = "RNA", layer = "data")
meta_data <- seurat_apical@meta.data
cellchat_apical <- createCellChat(object = data.input, meta = meta_data, group.by = "cellchat_labels")
# 5. データベースのセットアップと解析パイプラインの実行
CellChatDB <- CellChatDB.human
cellchat_apical@DB <- CellChatDB
cellchat_apical <- subsetData(cellchat_apical)
cellchat_apical <- identifyOverExpressedGenes(cellchat_apical)
# 1. どの列にちゃんとした細胞名（文字）が入っているか確認するため、主要な列の先頭数件を表示
head(seurat_apical@meta.data$cell_type)
head(seurat_apical@meta.data$seurat_clusters)
# 2. もし 'cell_type' にちゃんとした細胞種が入っている場合は、以下のようにグループを再指定して実行します
cellchat_apical <- createCellChat(object = data.input, meta = meta_data, group.by = "cell_type")
# 3. データベースの設定と再開
cellchat_apical@DB <- CellChatDB.human
cellchat_apical <- subsetData(cellchat_apical)
cellchat_apical <- identifyOverExpressedGenes(cellchat_apical)
cellchat_apical <- identifyOverExpressedInteractions(cellchat_apical)
# ここでエラーが出なければ成功です！続けて計算を進めます
cellchat_apical <- computeCommunProb(cellchat_apical, type = "triMean")
cellchat_apical <- filterCommunication(cellchat_apical, min.cells = 10)
cellchat_apical <- computeCommunProbPathway(cellchat_apical)
cellchat_apical <- aggregateNet(cellchat_apical)
# Figure 7 のバブルプロット出力
target_pathways <- c("COLLAGEN", "FN1", "LAMININ", "PTN", "APP", "MIF")
p_bubble <- netVisual_bubble(
cellchat_apical,
sources.use = "Macrophage",
targets.use = c("Fibroblast", "Epithelial"),
signaling = target_pathways,
remove.isolate = FALSE,
title.name = "Signaling from Macrophages in Apical Periodontitis"
) + theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
# 正しいシグナル経路のみを指定して再実行
target_pathways <- c("COLLAGEN", "FN1", "LAMININ", "APP", "MIF")
p_bubble <- netVisual_bubble(
cellchat_apical,
sources.use = "Macrophage",
targets.use = c("Fibroblast", "Epithelial"),
signaling = target_pathways,
remove.isolate = FALSE,
title.name = "Signaling from Macrophages in Apical Periodontitis"
) +
theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
# PDFとして保存
pdf("Figure7_CellChat_Macrophage_Signals.pdf", width = 7, height = 9)
print(p_bubble)
dev.off()
# 画面にも表示
print(p_bubble)
message("★ Figure 7 のバブルプロット出力が正常に完了しました！")
# オブジェクトに含まれている全ての細胞種を確認
levels(cellchat_apical@idents$joint)
# オブジェクトに含まれている全ての細胞種を確認
levels(cellchat_apical@idents$joint)
# 1. 念のため、Seuratオブジェクトの 'cell_type' 列にどのような細胞種が入っているか確認します
print(unique(seurat_apical$cell_type))
# 2. セルタイプをIdentsに再設定して、新しく綺麗なCellChatオブジェクトを作成します
Idents(seurat_apical) <- seurat_apical$cell_type
data.input <- GetAssayData(seurat_apical, assay = "RNA", layer = "data")
meta_data <- seurat_apical@meta.data
cellchat_apical <- createCellChat(object = data.input, meta = meta_data, group.by = "cell_type")
# 3. データベースのセットアップと計算パイプラインの一括実行
cellchat_apical@DB <- CellChatDB.human
cellchat_apical <- subsetData(cellchat_apical)
cellchat_apical <- identifyOverExpressedGenes(cellchat_apical)
cellchat_apical <- identifyOverExpressedInteractions(cellchat_apical)
cellchat_apical <- computeCommunProb(cellchat_apical, type = "triMean")
cellchat_apical <- filterCommunication(cellchat_apical, min.cells = 10)
cellchat_apical <- computeCommunProbPathway(cellchat_apical)
cellchat_apical <- aggregateNet(cellchat_apical)
# 4. オブジェクトに実際に含まれているパスウェイ名を確認して、確実に存在するものを指定します
print(cellchat_apical@netP$pathways)
# 5. Figure 7: バブルプロットの作成とPDF保存
# ※上の一覧に実際に含まれているパスウェイ名（例: "APP", "MIF", "LAMININ" 等）に合わせて指定してください
target_pathways <- intersect(c("COLLAGEN", "FN1", "LAMININ", "APP", "MIF"), cellchat_apical@netP$pathways)
p_bubble <- netVisual_bubble(
cellchat_apical,
sources.use = "Macrophage",
targets.use = c("Fibroblast", "Epithelial"),
signaling = target_pathways,
remove.isolate = FALSE,
title.name = "Signaling from Macrophages in Apical Periodontitis"
) +
theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
pdf("Figure7_CellChat_Macrophage_Signals.pdf", width = 7, height = 9)
print(p_bubble)
dev.off()
print(p_bubble)
message("★ Figure 7 の再描画が完了しました！")
# 統合された元のオブジェクト全体に含まれている細胞種を確認
print(unique(seurat_combined$cell_type))
# すでに統合・計算されている cellchat_merged に含まれる細胞種を確認する
levels(cellchat_merged@idents$joint)
# 1. cellchat_merged に含まれるクラスターの一覧を再確認（ご参考）
levels(cellchat_merged@idents$joint)
# 2. 利用可能なパスウェイの確認と共通部分の抽出
target_pathways <- intersect(c("COLLAGEN", "FN1", "LAMININ", "APP", "MIF"), cellchat_merged@netP$pathways)
# 3. クラスター番号を指定してバブルプロットを作成
# ※ "Cluster_3"（マクロファージ想定）から "Cluster_4", "Cluster_6"（線維芽細胞・上皮細胞想定）へのシグナル
# ※ 実際の細胞種に対応するクラスター番号に書き換えてご活用ください
p_bubble <- netVisual_bubble(
cellchat_merged,
sources.use = "Cluster_3",
targets.use = c("Cluster_4", "Cluster_6"),
signaling = target_pathways,
remove.isolate = FALSE,
title.name = "Signaling from Macrophage-related Cluster"
) +
theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
# cellchat_perio に含まれるクラスター・細胞グループの確認
levels(cellchat_perio@idents$joint)
# cellchat_perio に含まれるクラスター・細胞グループの確認
levels(cellchat_perio@idents$joint)
# 各メタデータ列に入っている値のサンプルを確認し、正しい細胞名・クラスター名がどこにあるか探します
sapply(seurat_combined@meta.data[, c("cell_type", "cellchat_labels", "seurat_clusters")], function(x) head(unique(x), 10))
# 1. 各クラスター（seurat_clusters）にどの細胞種が含まれているかを確認します
table(seurat_combined$seurat_clusters, seurat_combined$cell_type)
# 2. 確認した上で、適切なクラスター番号（例: マクロファージのクラスター と ターゲットのクラスター）を指定してCellChatオブジェクトを作成します
# ※以下は例です。上のテーブル結果を見て、マクロファージが含まれるクラスター番号（例: "3" など）に書き換えてください。
Idents(seurat_combined) <- seurat_combined$seurat_clusters
data.input <- GetAssayData(seurat_combined, assay = "RNA", layer = "data")
meta_data <- seurat_combined@meta.data
cellchat_complete <- createCellChat(object = data.input, meta = meta_data, group.by = "seurat_clusters")
# 3. パイプラインの実行
CellChatDB <- CellChatDB.human
cellchat_complete@DB <- CellChatDB
cellchat_complete <- subsetData(cellchat_complete)
cellchat_complete <- identifyOverExpressedGenes(cellchat_complete)
cellchat_complete <- identifyOverExpressedInteractions(cellchat_complete)
cellchat_complete <- computeCommunProb(cellchat_complete, type = "triMean")
cellchat_complete <- filterCommunication(cellchat_complete, min.cells = 10)
cellchat_complete <- computeCommunProbPathway(cellchat_complete)
cellchat_complete <- aggregateNet(cellchat_complete)
# 4. Figure 7: バブルプロットの作成とPDF保存
target_pathways <- intersect(c("COLLAGEN", "FN1", "LAMININ", "APP", "MIF"), cellchat_complete@netP$pathways)
# ※ sources.use と targets.use には、テーブルで確認した実際のクラスター番号（例: "3", "4" など）を指定します
p_bubble <- netVisual_bubble(
cellchat_complete,
sources.use = "3", # 例：マクロファージのクラスター番号
targets.use = c("4", "6"), # 例：ターゲットのクラスター番号
signaling = target_pathways,
remove.isolate = FALSE,
title.name = "Signaling from Macrophage Cluster"
) +
theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
pdf("Figure7_CellChat_Macrophage_Signals.pdf", width = 7, height = 9)
print(p_bubble)
dev.off()
print(p_bubble)
message("★ Figure 7 の作成が完了しました！")
# 環境内にあるSeurat関連のオブジェクトや、すべてのオブジェクト名を確認する
ls(pattern = "seurat")
# 主要なオブジェクトに含まれる cell_type の種類を一括で確認する
cat("--- seurat_combined ---\n"); print(unique(seurat_combined$cell_type))
cat("--- seurat_perio ---\n"); print(unique(seurat_perio$cell_type))
# seurat_perio と seurat_sub のメタデータの列名を確認する
colnames(seurat_perio@meta.data)
colnames(seurat_sub@meta.data)
# seurat_perio と seurat_sub のメタデータの列名を確認する
colnames(seurat_perio@meta.data)
colnames(seurat_sub@meta.data)
# seurat_perio の cellchat_labels に入っている内容を確認する
unique(seurat_perio$cellchat_labels)
# 1. cellchat_labels を用いて seurat_perio から CellChat オブジェクトを正しく作成
Idents(seurat_perio) <- seurat_perio$cellchat_labels
data.input <- GetAssayData(seurat_perio, assay = "RNA", layer = "data")
meta_data <- seurat_perio@meta.data
cellchat_perio_obj <- createCellChat(object = data.input, meta = meta_data, group.by = "cellchat_labels")
# 2. データベースのセットアップと計算パイプラインの一括実行
CellChatDB <- CellChatDB.human
cellchat_perio_obj@DB <- CellChatDB
cellchat_perio_obj <- subsetData(cellchat_perio_obj)
cellchat_perio_obj <- identifyOverExpressedGenes(cellchat_perio_obj)
cellchat_perio_obj <- identifyOverExpressedInteractions(cellchat_perio_obj)
cellchat_perio_obj <- computeCommunProb(cellchat_perio_obj, type = "triMean")
cellchat_perio_obj <- filterCommunication(cellchat_perio_obj, min.cells = 10)
cellchat_perio_obj <- computeCommunProbPathway(cellchat_perio_obj)
cellchat_perio_obj <- aggregateNet(cellchat_perio_obj)
# 3. Figure 7: バブルプロットの作成とPDF保存
target_pathways <- intersect(c("COLLAGEN", "FN1", "LAMININ", "APP", "MIF"), cellchat_perio_obj@netP$pathways)
# ※ sources.use にマクロファージが含まれるクラスター（例: "Cluster_3"）、
#    targets.use にターゲットとなるクラスター（例: "Cluster_4", "Cluster_6" 等）を指定します
p_bubble <- netVisual_bubble(
cellchat_perio_obj,
sources.use = "Cluster_3",
targets.use = c("Cluster_4", "Cluster_6"),
signaling = target_pathways,
remove.isolate = FALSE,
title.name = "Signaling from Macrophage Cluster in Apical Periodontitis"
) +
theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
# PDFとして保存
pdf("Figure7_CellChat_Macrophage_Signals.pdf", width = 7, height = 9)
print(p_bubble)
dev.off()
# 画面に表示
print(p_bubble)
message("★ Figure 7 のバブルプロット出力が正常に完了しました！")
# 絞り込みを行わず、オブジェクト内に存在するすべての有意なコミュニケーションを俯瞰するバブルプロット
p_all <- netVisual_bubble(cellchat_perio_obj, remove.isolate = TRUE)
# 画面に表示
print(p_all)
# 1. 縦軸の文字（リガンド・レセプター対）のフォントサイズを調整し、重なりを防ぐ
p_all_adjusted <- netVisual_bubble(cellchat_perio_obj, remove.isolate = TRUE) +
theme(
axis.text.y = element_text(size = 7), # 縦軸の文字サイズを小さくして重なりを解消
axis.text.x = element_text(angle = 45, hjust = 1, size = 8) # 横軸のクラスター名も見やすく傾ける
)
# 2. 縦長の高解像度PDFとして保存（高さ16インチに設定して文字が潰れないようにする）
pdf("Figure7_CellChat_All_Signals_Fixed.pdf", width = 10, height = 16)
print(p_all_adjusted)
dev.off()
# 3. 画面にも表示
print(p_all_adjusted)
message("★ 縦軸の重なりを解消した高解像度PDFの出力が完了しました！")
# 1. 縦軸の文字サイズを極限まで小さくし、行間や余白を最適化
p_all_extreme <- netVisual_bubble(cellchat_perio_obj, remove.isolate = TRUE) +
theme(
axis.text.y = element_text(size = 3.5, lineheight = 0.8), # 縦軸の文字サイズをさらに縮小
axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
)
# 2. 縦幅を「30インチ」という超縦長に設定してPDF出力
pdf("Figure7_CellChat_All_Signals_UltraTall.pdf", width = 12, height = 30)
print(p_all_extreme)
dev.off()
# 3. 画面にも表示
print(p_all_extreme)
message("★ 超縦長の高解像度PDFの出力が完了しました！")
# 1. 各細胞種の代表的なマーカー遺伝子をチェック
markers_to_check <- c("CD68", "LYZ", "COL1A1", "DCN", "KRT14", "KRT19", "CD3D", "NKG7")
# 2. クラスターごとの発現状況をドットプロットで可視化
DotPlot(seurat_perio, features = markers_to_check, group.by = "cellchat_labels") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
# マクロファージおよび上皮・線維芽細胞のマーカーで再確認
markers_target <- c("CD68", "CSF1R", "COL1A1", "DCN", "KRT14", "KRT19")
DotPlot(seurat_perio, features = markers_target, group.by = "cellchat_labels") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 1. ソースをマクロファージ（Cluster_7）、ターゲットを線維芽細胞・上皮細胞（Cluster_5, Cluster_4）に限定してバブルプロットを作成
p_fig7 <- netVisual_bubble(
cellchat_perio_obj,
sources.use = "Cluster_7",
targets.use = c("Cluster_5", "Cluster_4"),
remove.isolate = FALSE,
title.name = "Signaling from Macrophages to Fibroblasts & Epithelial Cells"
) +
theme(
plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
axis.text.y = element_text(size = 10)
)
# 2. 論文用の高解像度PDFとして保存
pdf("Figure7_Macrophage_to_Fibroblast_Epithelial.pdf", width = 8, height = 10)
print(p_fig7)
dev.off()
# 3. 画面にも表示
print(p_fig7)
message("★ 完璧な Figure 7 のPDF出力が完了しました！")
# 1. ソースをマクロファージ（Cluster_7）、ターゲットを線維芽細胞・上皮細胞（Cluster_5, Cluster_4）に限定してバブルプロットを作成
p_fig7 <- netVisual_bubble(
cellchat_perio_obj,
sources.use = "Cluster_7",
targets.use = c("Cluster_5", "Cluster_4"),
remove.isolate = FALSE,
title.name = "Signaling from Macrophages to Fibroblasts & Epithelial Cells"
) +
theme(
plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
axis.text.y = element_text(size = 10)
)
# 2. 論文用の高解像度PDFとして保存
pdf("Figure7_Macrophage_to_Fibroblast_Epithelial.pdf", width = 8, height = 10)
print(p_fig7)
dev.off()
# 3. 画面にも表示
print(p_fig7)
message("★ 完璧な Figure 7 のPDF出力が完了しました！")
# 1. 下側の余白（margin）を広げ、横軸のラベルが切れないように調整
p_fig7_fixed <- p_fig7 +
theme(
axis.text.x = element_text(angle = 45, hjust = 1, size = 10, vjust = 1),
plot.margin = margin(t = 10, r = 10, b = 40, l = 10) # 下側（b）の余白を大きく確保
)
# 2. 論文用の高解像度PDFとして再保存
pdf("Figure7_Macrophage_to_Fibroblast_Epithelial.pdf", width = 8, height = 10)
print(p_fig7_fixed)
dev.off()
# 3. 画面にも表示
print(p_fig7_fixed)
message("★ 横軸の文字切れを解消した完璧な Figure 7 の出力が完了しました！")
# 1. cellchat_perio_obj のオブジェクト内で、クラスタラベルの表記を細胞種名入りに変更する
# (Cluster_7 -> Macrophage, Cluster_4 -> Epithelial, Cluster_5 -> Fibroblast)
cellchat_perio_obj@net$centr # 念のため元の構造を保持しつつプロット時のラベルをカスタマイズします
# 2. ラベル名を書き換えたプロットの作成
p_fig7_labeled <- netVisual_bubble(
cellchat_perio_obj,
sources.use = "Cluster_7",
targets.use = c("Cluster_5", "Cluster_4"),
remove.isolate = FALSE,
title.name = "Signaling from Macrophages to Fibroblasts & Epithelial Cells"
)
# 横軸のラベル（因子水準）を直感的な名称に置換
# ※元々のラベルの並び順に合わせて書き換えます
current_labels <- levels(p_fig7_labeled$data$interaction_name_2) # またはx軸のラベル構造を確認
# gpplotのスケールを操作してラベル名を置換します
p_fig7_labeled <- p_fig7_labeled +
scale_x_discrete(
labels = c(
"7 -> Cluster_4" = "Macrophage (Cl. 7) -> Epithelial (Cl. 4)",
"7 -> Cluster_5" = "Macrophage (Cl. 7) -> Fibroblast (Cl. 5)"
)
) +
theme(
plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold"),
axis.text.y = element_text(size = 10),
plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
)
# 3. 論文用の高解像度PDFとして上書き保存
pdf("Figure7_Macrophage_to_Fibroblast_Epithelial.pdf", width = 9, height = 10)
print(p_fig7_labeled)
dev.off()
# 4. 画面にも表示
print(p_fig7_labeled)
message("★ 横軸に細胞種名を明記した Figure 7 の更新が完了しました！")
# 1. cellchat_perio_obj のオブジェクト内で、クラスタラベルの表記を細胞種名入りに変更する
# (Cluster_7 -> Macrophage, Cluster_4 -> Epithelial, Cluster_5 -> Fibroblast)
cellchat_perio_obj@net$centr # 念のため元の構造を保持しつつプロット時のラベルをカスタマイズします
# 2. ラベル名を書き換えたプロットの作成
p_fig7_labeled <- netVisual_bubble(
cellchat_perio_obj,
sources.use = "Cluster_7",
targets.use = c("Cluster_5", "Cluster_4"),
remove.isolate = FALSE,
title.name = "Signaling from Macrophages to Fibroblasts & Epithelial Cells"
)
# 横軸のラベル（因子水準）を直感的な名称に置換
# ※元々のラベルの並び順に合わせて書き換えます
current_labels <- levels(p_fig7_labeled$data$interaction_name_2) # またはx軸のラベル構造を確認
# gpplotのスケールを操作してラベル名を置換します
p_fig7_labeled <- p_fig7_labeled +
scale_x_discrete(
labels = c(
"7 -> Cluster_4" = "Macrophage (Cl. 7) -> Epithelial (Cl. 4)",
"7 -> Cluster_5" = "Macrophage (Cl. 7) -> Fibroblast (Cl. 5)"
)
) +
theme(
plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold"),
axis.text.y = element_text(size = 10),
plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
)
# 3. 論文用の高解像度PDFとして上書き保存
pdf("Figure7_Macrophage_to_Fibroblast_Epithelial.pdf", width = 9, height = 10)
print(p_fig7_labeled)
dev.off()
# 4. 画面にも表示
print(p_fig7_labeled)
message("★ 横軸に細胞種名を明記した Figure 7 の更新が完了しました！")
# 1. 横軸のラベルを改行（\n）を用いて、クラスター番号の下に細胞種名が来るように置換
p_fig7_multiline <- p_fig7 +
scale_x_discrete(
labels = c(
"7 -> Cluster_4" = "Cluster_7\n(Macrophage) -> Cluster_4\n(Epithelial cells)",
"7 -> Cluster_5" = "Cluster_7\n(Macrophage) -> Cluster_5\n(Fibroblasts)"
)
) +
theme(
plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, face = "bold"), # 改行に合わせて水平配置に調整
axis.text.y = element_text(size = 10),
plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
)
# 2. 論文用の高解像度PDFとして上書き保存
pdf("Figure7_Macrophage_to_Fibroblast_Epithelial.pdf", width = 8, height = 10)
print(p_fig7_multiline)
dev.off()
# 3. 画面にも表示
print(p_fig7_multiline)
message("★ 横軸に細胞種名を2行で明記した Figure 7 の更新が完了しました！")
# 1. 横軸のラベルを改行（\n）を用いて、クラスター番号の下に細胞種名が来るように置換
p_fig7_multiline <- p_fig7 +
scale_x_discrete(
labels = c(
"7 -> Cluster_4" = "Cluster_7\n(Macrophage) -> Cluster_4\n(Epithelial cells)",
"7 -> Cluster_5" = "Cluster_7\n(Macrophage) -> Cluster_5\n(Fibroblasts)"
)
) +
theme(
plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, face = "bold"), # 改行に合わせて水平配置に調整
axis.text.y = element_text(size = 10),
plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
)
# 2. 論文用の高解像度PDFとして上書き保存
pdf("Figure7_Macrophage_to_Fibroblast_Epithelial.pdf", width = 8, height = 10)
print(p_fig7_multiline)
dev.off()
# 3. 画面にも表示
print(p_fig7_multiline)
message("★ 横軸に細胞種名を2行で明記した Figure 7 の更新が完了しました！")
# 1. 正しいラベル名を指定してX軸のテキストを細胞種名入りに置換する
p_fig7_multiline <- p_fig7 +
scale_x_discrete(
labels = c(
"Cluster_7 -> Cluster_4" = "Cluster_7 (Macrophage)\n-> Cluster_4 (Epithelial cells)",
"Cluster_7 -> Cluster_5" = "Cluster_7 (Macrophage)\n-> Cluster_5 (Fibroblasts)"
)
) +
theme(
plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, face = "bold", lineheight = 1.2),
axis.text.y = element_text(size = 10),
plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
)
# 2. 論文用の高解像度PDFとして上書き保存
pdf("Figure7_Macrophage_to_Fibroblast_Epithelial.pdf", width = 8, height = 10)
print(p_fig7_multiline)
dev.off()
# 3. 画面にも表示
print(p_fig7_multiline)
message("★ 横軸の細胞種名表記の更新が完了しました！")
savehistory("~/EndoPerio_scRNA/scRNA_analysis.Rhistory")

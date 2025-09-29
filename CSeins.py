import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Circle, Rectangle
import warnings
warnings.filterwarnings('ignore')

class BreastCancerGenomeAnalyzer:
    def __init__(self):
        self.colors = {'A': '#FF6B6B', 'T': '#4ECDC4', 'C': '#45B7D1', 'G': '#F9A602'}
        self.genes_brca = ['BRCA1', 'BRCA2', 'TP53', 'PTEN', 'PALB2', 'CHEK2', 'ATM', 'CDH1']
        
    def generate_genomic_data(self):
        """Génère des données génomiques simulées pour le cancer du sein"""
        print("🧬 Génération des données génomiques du cancer du sein...")
        
        # Mutations communes dans le cancer du sein
        mutations_data = {
            'Gene': self.genes_brca + ['PIK3CA', 'AKT1', 'GATA3', 'MAP3K1', 'ESR1', 'ERBB2'],
            'Mutation_Frequency': [0.25, 0.20, 0.35, 0.08, 0.12, 0.09, 0.07, 0.08, 
                                  0.28, 0.04, 0.11, 0.13, 0.06, 0.15],
            'Mutation_Type': ['Frameshift', 'Nonsense', 'Missense', 'Deletion', 
                            'Frameshift', 'Missense', 'Missense', 'Splice_Site',
                            'Missense', 'Missense', 'Missense', 'Frameshift', 
                            'Missense', 'Amplification'],
            'Clinical_Significance': ['Pathogenic', 'Pathogenic', 'Pathogenic', 'Pathogenic',
                                    'Pathogenic', 'Pathogenic', 'Pathogenic', 'Pathogenic',
                                    'Oncogenic', 'Oncogenic', 'Oncogenic', 'Oncogenic',
                                    'Oncogenic', 'Oncogenic']
        }
        
        return pd.DataFrame(mutations_data)
    
    def create_ctag_diagram(self, df):
        """Crée un diagramme CTAG pour visualiser les mutations"""
        plt.style.use('seaborn-v0_8')
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 16))
        
        # 1. Fréquence des mutations par gène
        self._plot_mutation_frequency(df, ax1)
        
        # 2. Répartition des types de mutations
        self._plot_mutation_types(df, ax2)
        
        # 3. Diagramme CTAG sequence
        self._plot_ctag_sequence(ax3)
        
        # 4. Impact clinique
        self._plot_clinical_impact(df, ax4)
        
        plt.suptitle('Analyse Génomique du Cancer du Sein - Diagramme CTAG', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig('breast_cancer_ctag_diagram.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Générer le rapport d'analyse
        self._generate_genomic_report(df)
    
    def _plot_mutation_frequency(self, df, ax):
        """Plot de la fréquence des mutations par gène"""
        colors = ['#FF6B6B' if 'BRCA' in gene else '#45B7D1' for gene in df['Gene']]
        
        bars = ax.barh(df['Gene'], df['Mutation_Frequency']*100, 
                      color=colors, alpha=0.8)
        
        # Ajouter les valeurs sur les barres
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 1, bar.get_y() + bar.get_height()/2, 
                   f'{width:.1f}%', ha='left', va='center', fontweight='bold')
        
        ax.set_xlabel('Fréquence des Mutations (%)')
        ax.set_title('Fréquence des Mutations par Gène\n dans le Cancer du Sein', 
                    fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='x')
    
    def _plot_mutation_types(self, df, ax):
        """Plot de la répartition des types de mutations"""
        mutation_counts = df['Mutation_Type'].value_counts()
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#F9A602', '#6A0572', '#2A9D8F']
        wedges, texts, autotexts = ax.pie(mutation_counts.values, 
                                         labels=mutation_counts.index,
                                         colors=colors[:len(mutation_counts)],
                                         autopct='%1.1f%%',
                                         startangle=90)
        
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
        
        ax.set_title('Répartition des Types de Mutations', 
                    fontsize=12, fontweight='bold')
    
    def _plot_ctag_sequence(self, ax):
        """Crée un diagramme CTAG sequence"""
        # Séquence simulée avec mutations
        sequence = "ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
        mutations = [5, 12, 25, 38, 45]  # positions de mutations
        
        ax.set_xlim(0, len(sequence))
        ax.set_ylim(0, 2)
        
        # Dessiner la séquence de base
        for i, base in enumerate(sequence):
            color = self.colors[base]
            rect = Rectangle((i, 0.8), 1, 0.4, facecolor=color, alpha=0.7, edgecolor='black')
            ax.add_patch(rect)
            ax.text(i + 0.5, 1, base, ha='center', va='center', 
                   fontweight='bold', fontsize=8)
        
        # Marquer les mutations
        for mut_pos in mutations:
            circle = Circle((mut_pos + 0.5, 1.6), 0.3, 
                          facecolor='red', alpha=0.8, edgecolor='darkred')
            ax.add_patch(circle)
            ax.text(mut_pos + 0.5, 1.6, 'M', ha='center', va='center', 
                   fontweight='bold', color='white', fontsize=8)
        
        ax.set_title('Diagramme CTAG - Séquence Génomique avec Mutations', 
                    fontsize=12, fontweight='bold')
        ax.set_xlabel('Position dans la séquence')
        ax.set_ylabel('')
        ax.set_yticks([])
        ax.grid(True, alpha=0.3)
    
    def _plot_clinical_impact(self, df, ax):
        """Plot de l'impact clinique des mutations"""
        clinical_counts = df['Clinical_Significance'].value_counts()
        
        colors = ['#E76F51', '#2A9D8F', '#F9A602']
        bars = ax.bar(clinical_counts.index, clinical_counts.values, 
                     color=colors[:len(clinical_counts)], alpha=0.8)
        
        # Ajouter les valeurs sur les barres
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.1,
                   f'{height}', ha='center', va='bottom', fontweight='bold')
        
        ax.set_ylabel('Nombre de Gènes')
        ax.set_title('Impact Clinique des Mutations', 
                    fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Rotation des labels x
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    
    def create_advanced_genomic_analysis(self, df):
        """Crée une analyse génomique avancée"""
        fig = plt.figure(figsize=(18, 12))
        
        # 1. Heatmap des mutations
        ax1 = plt.subplot(2, 3, 1)
        self._plot_mutation_heatmap(df, ax1)
        
        # 2. Distribution des fréquences
        ax2 = plt.subplot(2, 3, 2)
        self._plot_frequency_distribution(df, ax2)
        
        # 3. Réseau d'interaction des gènes
        ax3 = plt.subplot(2, 3, 3)
        self._plot_gene_interaction_network(ax3)
        
        # 4. Profil mutationnel
        ax4 = plt.subplot(2, 3, 4)
        self._plot_mutational_signature(ax4)
        
        # 5. Pathway analysis
        ax5 = plt.subplot(2, 3, 5)
        self._plot_pathway_analysis(df, ax5)
        
        # 6. CTAG detailed
        ax6 = plt.subplot(2, 3, 6)
        self._plot_detailed_ctag(ax6)
        
        plt.suptitle('Analyse Génomique Avancée - Cancer du Sein', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig('advanced_breast_cancer_genomics.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def _plot_mutation_heatmap(self, df, ax):
        """Crée une heatmap des mutations"""
        # Créer une matrice simulée pour la heatmap
        genes = df['Gene'].tolist()
        mutation_types = df['Mutation_Type'].unique()
        
        # Matrice de fréquence simulée
        heatmap_data = np.random.rand(len(genes), len(mutation_types)) * 0.3
        
        # Renforcer certaines associations
        for i, gene in enumerate(genes):
            if 'BRCA' in gene:
                heatmap_data[i, 0] = 0.8  # Forte association avec frameshift
            if 'PIK3' in gene:
                heatmap_data[i, 2] = 0.9  # Forte association avec missense
        
        sns.heatmap(heatmap_data, ax=ax, cmap='YlOrRd', 
                   xticklabels=mutation_types, yticklabels=genes,
                   cbar_kws={'label': 'Fréquence relative'})
        
        ax.set_title('Heatmap des Mutations\npar Gène et Type', fontsize=10, fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.tick_params(axis='y', rotation=0)
    
    def _plot_frequency_distribution(self, df, ax):
        """Plot de la distribution des fréquences de mutations"""
        ax.hist(df['Mutation_Frequency']*100, bins=10, 
               color='#45B7D1', alpha=0.7, edgecolor='black')
        
        ax.axvline(df['Mutation_Frequency'].mean()*100, color='red', 
                  linestyle='--', linewidth=2, label=f'Moyenne: {df["Mutation_Frequency"].mean()*100:.1f}%')
        
        ax.set_xlabel('Fréquence des Mutations (%)')
        ax.set_ylabel('Nombre de Gènes')
        ax.set_title('Distribution des Fréquences de Mutations', 
                    fontsize=10, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    def _plot_gene_interaction_network(self, ax):
        """Crée un réseau d'interaction des gènes"""
        # Positions simulées pour le réseau
        nodes = {
            'BRCA1': (2, 3), 'BRCA2': (4, 3), 'TP53': (3, 1),
            'PTEN': (1, 2), 'PALB2': (5, 2), 'CHEK2': (2, 4),
            'ATM': (4, 4), 'PIK3CA': (3, 5)
        }
        
        # Connexions simulées
        edges = [('BRCA1', 'BRCA2'), ('BRCA1', 'TP53'), ('BRCA2', 'PALB2'),
                ('TP53', 'PTEN'), ('BRCA1', 'CHEK2'), ('BRCA2', 'ATM'),
                ('PIK3CA', 'TP53'), ('CHEK2', 'ATM')]
        
        # Dessiner les arêtes
        for edge in edges:
            start = nodes[edge[0]]
            end = nodes[edge[1]]
            ax.plot([start[0], end[0]], [start[1], end[1]], 
                   'gray', alpha=0.6, linewidth=2)
        
        # Dessiner les nœuds
        for gene, pos in nodes.items():
            color = '#FF6B6B' if 'BRCA' in gene else '#45B7D1'
            ax.scatter(pos[0], pos[1], s=300, c=color, alpha=0.8, 
                      edgecolors='black', linewidth=2)
            ax.text(pos[0], pos[1], gene, ha='center', va='center', 
                   fontweight='bold', fontsize=8)
        
        ax.set_xlim(0, 6)
        ax.set_ylim(0, 6)
        ax.set_title('Réseau d\'Interaction des Gènes', 
                    fontsize=10, fontweight='bold')
        ax.set_xticks([])
        ax.set_yticks([])
    
    def _plot_mutational_signature(self, ax):
        """Plot des signatures mutationnelles"""
        # Signatures mutationnelles simulées
        contexts = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
        signature1 = [0.1, 0.05, 0.4, 0.1, 0.2, 0.15]  # Signature liée à l'âge
        signature2 = [0.3, 0.1, 0.2, 0.2, 0.1, 0.1]   # Signature BRCA
        
        x = np.arange(len(contexts))
        width = 0.35
        
        ax.bar(x - width/2, signature1, width, label='Signature Âge', 
               color='#4ECDC4', alpha=0.8)
        ax.bar(x + width/2, signature2, width, label='Signature BRCA', 
               color='#FF6B6B', alpha=0.8)
        
        ax.set_xlabel('Type de Mutation')
        ax.set_ylabel('Fréquence relative')
        ax.set_title('Signatures Mutationnelles', fontsize=10, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(contexts)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    def _plot_pathway_analysis(self, df, ax):
        """Analyse des voies de signalisation"""
        pathways = {
            'Réparation ADN': ['BRCA1', 'BRCA2', 'PALB2', 'ATM', 'CHEK2'],
            'PI3K/AKT': ['PIK3CA', 'PTEN', 'AKT1'],
            'Cycle Cellulaire': ['TP53', 'CDH1'],
            'Signalisation Hormonale': ['ESR1', 'GATA3']
        }
        
        pathway_counts = {pathway: len(genes) for pathway, genes in pathways.items()}
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#F9A602']
        bars = ax.bar(pathway_counts.keys(), pathway_counts.values(), 
                     color=colors, alpha=0.8)
        
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.1,
                   f'{height}', ha='center', va='bottom', fontweight='bold')
        
        ax.set_ylabel('Nombre de Gènes')
        ax.set_title('Analyse des Voies de Signalisation', 
                    fontsize=10, fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, axis='y')
    
    def _plot_detailed_ctag(self, ax):
        """Diagramme CTAG détaillé avec transitions/transversions"""
        # Données simulées pour les transitions/transversions
        mutation_types = ['C>T', 'C>G', 'C>A', 'T>C', 'T>G', 'T>A']
        counts = [150, 45, 60, 80, 35, 55]  # Comptes simulés
        
        colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#F9A602', '#6A0572', '#2A9D8F']
        bars = ax.bar(mutation_types, counts, color=colors, alpha=0.8)
        
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 2,
                   f'{height}', ha='center', va='bottom', fontweight='bold')
        
        ax.set_ylabel('Nombre de Mutations')
        ax.set_title('Spectre des Mutations CTAG\n(Transitions/Transversions)', 
                    fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
    
    def _generate_genomic_report(self, df):
        """Génère un rapport d'analyse génomique"""
        print("\n🧬 RAPPORT D'ANALYSE GÉNOMIQUE - CANCER DU SEIN")
        print("=" * 60)
        
        # 1. Gènes les plus fréquemment mutés
        print("\n1. 📊 GÈNES LES PLUS MUTÉS:")
        top_genes = df.nlargest(5, 'Mutation_Frequency')
        for _, row in top_genes.iterrows():
            print(f"   • {row['Gene']}: {row['Mutation_Frequency']*100:.1f}% "
                  f"({row['Mutation_Type']}) - {row['Clinical_Significance']}")
        
        # 2. Analyse BRCA
        print("\n2. 🔍 ANALYSE DES GÈNES BRCA:")
        brca_genes = df[df['Gene'].str.contains('BRCA')]
        for _, row in brca_genes.iterrows():
            print(f"   • {row['Gene']}: Fréquence {row['Mutation_Frequency']*100:.1f}%")
        
        # 3. Implications thérapeutiques
        print("\n3. 💊 IMPLICATIONS THÉRAPEUTIQUES:")
        print("   • Inhibiteurs de PARP: Indiqués pour mutations BRCA1/2")
        print("   • Inhibiteurs de PI3K: Pour mutations PIK3CA")
        print("   • Immunothérapie: Charge mutationnelle élevée")
        print("   • Thérapies ciblées: Basées sur le profil mutationnel")
        
        # 4. Recommandations de dépistage
        print("\n4. 🎯 RECOMMANDATIONS DE DÉPISTAGE:")
        print("   • Séquençage panel NGS: BRCA1, BRCA2, TP53, PALB2, etc.")
        print("   • Conseil génétique: Pour antécédents familiaux")
        print("   • Surveillance renforcée: Porteurs mutations à haut risque")
        
        # 5. Statistiques globales
        print("\n5. 📈 STATISTIQUES GLOBALES:")
        print(f"   • Nombre total de gènes analysés: {len(df)}")
        print(f"   • Fréquence moyenne des mutations: {df['Mutation_Frequency'].mean()*100:.1f}%")
        print(f"   • Gènes pathogènes: {len(df[df['Clinical_Significance'] == 'Pathogenic'])}")
        print(f"   • Gènes oncogéniques: {len(df[df['Clinical_Significance'] == 'Oncogenic'])}")

def main():
    """Fonction principale"""
    print("🧬 ANALYSE GÉNOMIQUE DU CANCER DU SEIN - DIAGRAMME CTAG")
    print("=" * 60)
    
    # Initialiser l'analyseur
    analyzer = BreastCancerGenomeAnalyzer()
    
    # Générer les données
    genomic_data = analyzer.generate_genomic_data()
    
    # Aperçu des données
    print("\n👀 Aperçu des données génomiques:")
    print(genomic_data.head(8))
    
    # Créer le diagramme CTAG principal
    print("\n📊 Création du diagramme CTAG...")
    analyzer.create_ctag_diagram(genomic_data)
    
    # Créer l'analyse avancée
    print("\n🔬 Création de l'analyse génomique avancée...")
    analyzer.create_advanced_genomic_analysis(genomic_data)
    
    print("\n✅ Analyse génomique terminée!")
    print("📁 Fichiers générés:")
    print("   • breast_cancer_ctag_diagram.png")
    print("   • advanced_breast_cancer_genomics.png")

if __name__ == "__main__":
    main()
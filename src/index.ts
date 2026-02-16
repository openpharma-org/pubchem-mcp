#!/usr/bin/env node
import { Server } from '@modelcontextprotocol/sdk/server/index.js';
import { StdioServerTransport } from '@modelcontextprotocol/sdk/server/stdio.js';
import {
  CallToolRequestSchema,
  ErrorCode,
  ListResourcesRequestSchema,
  ListResourceTemplatesRequestSchema,
  ListToolsRequestSchema,
  McpError,
  ReadResourceRequestSchema,
} from '@modelcontextprotocol/sdk/types.js';
import axios, { AxiosInstance } from 'axios';

// PubChem API interfaces
interface CompoundSearchResult {
  IdentifierList: {
    CID: number[];
  };
}

interface PropertyData {
  PropertyTable: {
    Properties: Array<{
      CID: number;
      MolecularFormula?: string;
      MolecularWeight?: number;
      CanonicalSMILES?: string;
      IsomericSMILES?: string;
      InChI?: string;
      InChIKey?: string;
      IUPACName?: string;
      XLogP?: number;
      TPSA?: number;
      Complexity?: number;
      Charge?: number;
      HBondDonorCount?: number;
      HBondAcceptorCount?: number;
      RotatableBondCount?: number;
      HeavyAtomCount?: number;
      AtomStereoCount?: number;
      DefinedAtomStereoCount?: number;
      BondStereoCount?: number;
      DefinedBondStereoCount?: number;
      Volume3D?: number;
      ConformerCount3D?: number;
    }>;
  };
}

// Type guards and validation functions
const isValidCompoundSearchArgs = (
  args: any
): args is { query: string; search_type?: string; max_records?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.query === 'string' &&
    args.query.length > 0 &&
    (args.search_type === undefined || ['name', 'smiles', 'inchi', 'sdf', 'cid', 'formula'].includes(args.search_type)) &&
    (args.max_records === undefined || (typeof args.max_records === 'number' && args.max_records > 0 && args.max_records <= 10000))
  );
};

const isValidCidArgs = (
  args: any
): args is { cid: number | string; format?: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    (typeof args.cid === 'number' || typeof args.cid === 'string') &&
    (args.format === undefined || ['json', 'sdf', 'xml', 'asnt', 'asnb'].includes(args.format))
  );
};

const isValidSmilesArgs = (
  args: any
): args is { smiles: string; threshold?: number; max_records?: number } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    typeof args.smiles === 'string' &&
    args.smiles.length > 0 &&
    (args.threshold === undefined || (typeof args.threshold === 'number' && args.threshold >= 0 && args.threshold <= 100)) &&
    (args.max_records === undefined || (typeof args.max_records === 'number' && args.max_records > 0 && args.max_records <= 10000))
  );
};

const isValidBatchArgs = (
  args: any
): args is { cids: number[]; operation?: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    Array.isArray(args.cids) &&
    args.cids.length > 0 &&
    args.cids.length <= 200 &&
    args.cids.every((cid: any) => typeof cid === 'number' && cid > 0) &&
    (args.operation === undefined || ['property', 'synonyms', 'classification', 'description'].includes(args.operation))
  );
};

const isValidConformerArgs = (
  args: any
): args is { cid: number | string; conformer_type?: string } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    (typeof args.cid === 'number' || typeof args.cid === 'string') &&
    (args.conformer_type === undefined || ['3d', '2d'].includes(args.conformer_type))
  );
};

const isValidPropertiesArgs = (
  args: any
): args is { cid: number | string; properties?: string[] } => {
  return (
    typeof args === 'object' &&
    args !== null &&
    (typeof args.cid === 'number' || typeof args.cid === 'string') &&
    (args.properties === undefined || (Array.isArray(args.properties) && args.properties.every((p: any) => typeof p === 'string')))
  );
};

class PubChemServer {
  private server: Server;
  private apiClient: AxiosInstance;

  constructor() {
    this.server = new Server(
      {
        name: 'pubchem-server',
        version: '1.0.0',
      },
      {
        capabilities: {
          resources: {},
          tools: {},
        },
      }
    );

    // Initialize PubChem API client
    this.apiClient = axios.create({
      baseURL: 'https://pubchem.ncbi.nlm.nih.gov/rest/pug',
      timeout: 30000,
      headers: {
        'User-Agent': 'PubChem-MCP-Server/1.0.0',
        'Accept': 'application/json',
      },
    });

    this.setupResourceHandlers();
    this.setupToolHandlers();

    // Error handling
    this.server.onerror = (error: any) => console.error('[MCP Error]', error);
    process.on('SIGINT', async () => {
      await this.server.close();
      process.exit(0);
    });
  }

  private setupResourceHandlers() {
    // List available resource templates
    this.server.setRequestHandler(
      ListResourceTemplatesRequestSchema,
      async () => ({
        resourceTemplates: [
          {
            uriTemplate: 'pubchem://compound/{cid}',
            name: 'PubChem compound entry',
            mimeType: 'application/json',
            description: 'Complete compound information for a PubChem CID',
          },
          {
            uriTemplate: 'pubchem://structure/{cid}',
            name: 'Chemical structure data',
            mimeType: 'application/json',
            description: '2D/3D structure information for a compound',
          },
          {
            uriTemplate: 'pubchem://properties/{cid}',
            name: 'Chemical properties',
            mimeType: 'application/json',
            description: 'Molecular properties and descriptors for a compound',
          },
          {
            uriTemplate: 'pubchem://bioassay/{aid}',
            name: 'PubChem bioassay data',
            mimeType: 'application/json',
            description: 'Bioassay information and results for an AID',
          },
          {
            uriTemplate: 'pubchem://similarity/{smiles}',
            name: 'Similarity search results',
            mimeType: 'application/json',
            description: 'Chemical similarity search results for a SMILES string',
          },
          {
            uriTemplate: 'pubchem://safety/{cid}',
            name: 'Safety and toxicity data',
            mimeType: 'application/json',
            description: 'Safety classifications and toxicity information',
          },
        ],
      })
    );

    // Handle resource requests
    this.server.setRequestHandler(
      ReadResourceRequestSchema,
      async (request: any) => {
        const uri = request.params.uri;

        // Handle compound info requests
        const compoundMatch = uri.match(/^pubchem:\/\/compound\/([0-9]+)$/);
        if (compoundMatch) {
          const cid = compoundMatch[1];
          try {
            const response = await this.apiClient.get(`/compound/cid/${cid}/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch compound ${cid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle structure requests
        const structureMatch = uri.match(/^pubchem:\/\/structure\/([0-9]+)$/);
        if (structureMatch) {
          const cid = structureMatch[1];
          try {
            const response = await this.apiClient.get(`/compound/cid/${cid}/property/CanonicalSMILES,IsomericSMILES,InChI,InChIKey/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch structure for ${cid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle properties requests
        const propertiesMatch = uri.match(/^pubchem:\/\/properties\/([0-9]+)$/);
        if (propertiesMatch) {
          const cid = propertiesMatch[1];
          try {
            const response = await this.apiClient.get(`/compound/cid/${cid}/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,Complexity/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch properties for ${cid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle bioassay requests
        const bioassayMatch = uri.match(/^pubchem:\/\/bioassay\/([0-9]+)$/);
        if (bioassayMatch) {
          const aid = bioassayMatch[1];
          try {
            const response = await this.apiClient.get(`/assay/aid/${aid}/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch bioassay ${aid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle similarity search requests
        const similarityMatch = uri.match(/^pubchem:\/\/similarity\/(.+)$/);
        if (similarityMatch) {
          const smiles = decodeURIComponent(similarityMatch[1]);
          try {
            const response = await this.apiClient.post('/compound/similarity/smiles/JSON', {
              smiles: smiles,
              Threshold: 90,
              MaxRecords: 100,
            });
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to perform similarity search: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        // Handle safety data requests
        const safetyMatch = uri.match(/^pubchem:\/\/safety\/([0-9]+)$/);
        if (safetyMatch) {
          const cid = safetyMatch[1];
          try {
            const response = await this.apiClient.get(`/compound/cid/${cid}/classification/JSON`);
            return {
              contents: [
                {
                  uri: request.params.uri,
                  mimeType: 'application/json',
                  text: JSON.stringify(response.data, null, 2),
                },
              ],
            };
          } catch (error) {
            throw new McpError(
              ErrorCode.InternalError,
              `Failed to fetch safety data for ${cid}: ${error instanceof Error ? error.message : 'Unknown error'}`
            );
          }
        }

        throw new McpError(
          ErrorCode.InvalidRequest,
          `Invalid URI format: ${uri}`
        );
      }
    );
  }

  private setupToolHandlers() {
    this.server.setRequestHandler(ListToolsRequestSchema, async () => ({
      tools: [
        {
          name: 'pubchem',
          description: 'Unified PubChem database access tool for chemical compound information, structure analysis, molecular properties, bioassay data, and safety information. Provides access to over 110 million compounds with extensive chemical informatics capabilities.',
          inputSchema: {
            type: 'object',
            properties: {
              method: {
                type: 'string',
                enum: [
                  'search_compounds',
                  'get_compound_info',
                  'search_by_smiles',
                  'get_compound_synonyms',
                  'search_similar_compounds',
                  'get_3d_conformers',
                  'analyze_stereochemistry',
                  'get_compound_properties',
                  'get_assay_info',
                  'get_safety_data',
                  'batch_compound_lookup',
                  'get_patent_ids',
                  // Unimplemented methods (commented out):
                  // 'search_by_inchi',
                  // 'search_by_cas_number',
                  // 'substructure_search',
                  // 'superstructure_search',
                  // 'calculate_descriptors',
                  // 'predict_admet_properties',
                  // 'assess_drug_likeness',
                  // 'analyze_molecular_complexity',
                  // 'get_pharmacophore_features',
                  // 'search_bioassays',
                  // 'get_compound_bioactivities',
                  // 'search_by_target',
                  // 'compare_activity_profiles',
                  // 'get_toxicity_info',
                  // 'assess_environmental_fate',
                  // 'get_regulatory_info',
                  // 'get_external_references',
                  // 'search_patents',
                  // 'get_literature_references',
                ],
                description: 'The PubChem operation to perform: search_compounds (search by name/CAS/formula), get_compound_info (detailed compound by CID), search_by_smiles (exact SMILES match), get_compound_synonyms (all names), search_similar_compounds (Tanimoto similarity), get_3d_conformers (3D structural data), analyze_stereochemistry (chirality analysis), get_compound_properties (MW/logP/TPSA), get_assay_info (detailed assay by AID), get_safety_data (GHS classifications), batch_compound_lookup (bulk processing up to 200 compounds), get_patent_ids (patent IDs associated with a compound by CID or SMILES)',
              },
              query: {
                type: 'string',
                description: 'Search query string (for search methods)',
              },
              cid: {
                type: ['number', 'string'],
                description: 'PubChem Compound ID (CID)',
              },
              aid: {
                type: 'number',
                description: 'PubChem Assay ID (AID)',
              },
              smiles: {
                type: 'string',
                description: 'SMILES string of the query molecule',
              },
              inchi: {
                type: 'string',
                description: 'InChI string or InChI key',
              },
              cas_number: {
                type: 'string',
                description: 'CAS Registry Number (e.g., 50-78-2)',
              },
              search_type: {
                type: 'string',
                enum: ['name', 'smiles', 'inchi', 'sdf', 'cid', 'formula'],
                description: 'Type of search to perform (default: name)',
              },
              max_records: {
                type: 'number',
                description: 'Maximum number of results (1-10000, default: 100)',
                minimum: 1,
                maximum: 10000,
              },
              threshold: {
                type: 'number',
                description: 'Similarity threshold (0-100, default: 90)',
                minimum: 0,
                maximum: 100,
              },
              properties: {
                type: 'array',
                items: { type: 'string' },
                description: 'Specific properties to retrieve',
              },
              format: {
                type: 'string',
                enum: ['json', 'sdf', 'xml', 'asnt', 'asnb'],
                description: 'Output format (default: json)',
              },
              conformer_type: {
                type: 'string',
                enum: ['3d', '2d'],
                description: 'Type of conformer data (default: 3d)',
              },
              descriptor_type: {
                type: 'string',
                enum: ['all', 'basic', 'topological', '3d'],
                description: 'Type of descriptors (default: all)',
              },
              target: {
                type: 'string',
                description: 'Target name (gene, protein, or pathway)',
              },
              activity_type: {
                type: 'string',
                description: 'Type of activity (e.g., IC50, EC50, Ki)',
              },
              activity_outcome: {
                type: 'string',
                enum: ['active', 'inactive', 'inconclusive', 'all'],
                description: 'Filter by activity outcome (default: all)',
              },
              cids: {
                type: 'array',
                items: { type: 'number' },
                description: 'Array of PubChem CIDs (for batch/comparison operations)',
              },
              operation: {
                type: 'string',
                enum: ['property', 'synonyms', 'classification', 'description'],
                description: 'Batch operation to perform (default: property)',
              },
              source: {
                type: 'string',
                description: 'Data source (e.g., ChEMBL, NCGC)',
              },
            },
            required: ['method'],
          },
        },
      ],
    }));

    this.server.setRequestHandler(CallToolRequestSchema, async (request: any) => {
      const { name, arguments: args } = request.params;

      if (name !== 'pubchem') {
        throw new McpError(
          ErrorCode.MethodNotFound,
          `Unknown tool: ${name}`
        );
      }

      if (!args.method || typeof args.method !== 'string') {
        throw new McpError(
          ErrorCode.InvalidParams,
          'The "method" parameter is required and must be a string'
        );
      }

      try {
        let result;
        switch (args.method) {
          // Chemical Search & Retrieval
          case 'search_compounds':
            result = await this.handleSearchCompounds(args);
            break;
          case 'get_compound_info':
            result = await this.handleGetCompoundInfo(args);
            break;
          case 'search_by_smiles':
            result = await this.handleSearchBySmiles(args);
            break;
          case 'get_compound_synonyms':
            result = await this.handleGetCompoundSynonyms(args);
            break;

          // Structure Analysis & Similarity
          case 'search_similar_compounds':
            result = await this.handleSearchSimilarCompounds(args);
            break;
          case 'get_3d_conformers':
            result = await this.handleGet3dConformers(args);
            break;
          case 'analyze_stereochemistry':
            result = await this.handleAnalyzeStereochemistry(args);
            break;

          // Chemical Properties & Descriptors
          case 'get_compound_properties':
            result = await this.handleGetCompoundProperties(args);
            break;

          // Bioassay & Activity Data
          case 'get_assay_info':
            result = await this.handleGetAssayInfo(args);
            break;

          // Safety & Toxicity
          case 'get_safety_data':
            result = await this.handleGetSafetyData(args);
            break;

          // Cross-References & Integration
          case 'batch_compound_lookup':
            result = await this.handleBatchCompoundLookup(args);
            break;

          // Patents
          case 'get_patent_ids':
            result = await this.handleGetPatentIds(args);
            break;

          default:
            throw new McpError(
              ErrorCode.MethodNotFound,
              `Unknown method: ${args.method}`
            );
        }

        return result;
      } catch (error) {
        return {
          content: [
            {
              type: 'text',
              text: `Error executing method ${args.method}: ${error instanceof Error ? error.message : 'Unknown error'}`,
            },
          ],
          isError: true,
        };
      }
    });
  }

  // Chemical Search & Retrieval handlers
  private async handleSearchCompounds(args: any) {
    if (!isValidCompoundSearchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid compound search arguments');
    }

    try {
      const searchType = args.search_type || 'name';
      const maxRecords = args.max_records || 100;

      const response = await this.apiClient.get(`/compound/${searchType}/${encodeURIComponent(args.query)}/cids/JSON`, {
        params: {
          MaxRecords: maxRecords,
        },
      });

      if (response.data?.IdentifierList?.CID?.length > 0) {
        const cids = response.data.IdentifierList.CID.slice(0, 10);
        const detailsResponse = await this.apiClient.get(`/compound/cid/${cids.join(',')}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON`);

        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                query: args.query,
                search_type: searchType,
                total_found: response.data.IdentifierList.CID.length,
                details: detailsResponse.data,
              }, null, 2),
            },
          ],
        };
      }

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({ message: 'No compounds found', query: args.query }, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search compounds: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleGetCompoundInfo(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid CID arguments');
    }

    try {
      const format = args.format || 'json';
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/${format === 'json' ? 'JSON' : format}`);

      return {
        content: [
          {
            type: 'text',
            text: format === 'json' ? JSON.stringify(response.data, null, 2) : String(response.data),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get compound info: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleSearchBySmiles(args: any) {
    if (!isValidSmilesArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid SMILES arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/smiles/${encodeURIComponent(args.smiles)}/cids/JSON`);

      if (response.data?.IdentifierList?.CID?.length > 0) {
        const cid = response.data.IdentifierList.CID[0];
        const detailsResponse = await this.apiClient.get(`/compound/cid/${cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON`);

        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify({
                query_smiles: args.smiles,
                found_cid: cid,
                details: detailsResponse.data,
              }, null, 2),
            },
          ],
        };
      }

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({ message: 'No exact match found', query_smiles: args.smiles }, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search by SMILES: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  // Placeholder implementations (commented out - uncomment to implement)
  // private async handleSearchByInchi(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'InChI search not yet implemented', args }, null, 2) }] };
  // }

  // private async handleSearchByCasNumber(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'CAS search not yet implemented', args }, null, 2) }] };
  // }

  private async handleGetCompoundSynonyms(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid CID arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/synonyms/JSON`);
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get compound synonyms: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleSearchSimilarCompounds(args: any) {
    if (!isValidSmilesArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid similarity search arguments');
    }

    try {
      const threshold = args.threshold || 90;
      const maxRecords = args.max_records || 100;

      const response = await this.apiClient.post('/compound/similarity/smiles/JSON', {
        smiles: args.smiles,
        Threshold: threshold,
        MaxRecords: maxRecords,
      });

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to search similar compounds: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  // private async handleSubstructureSearch(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Substructure search not yet implemented', args }, null, 2) }] };
  // }

  // private async handleSuperstructureSearch(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Superstructure search not yet implemented', args }, null, 2) }] };
  // }

  private async handleGet3dConformers(args: any) {
    if (!isValidConformerArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid 3D conformer arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/property/Volume3D,ConformerCount3D/JSON`);

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({
              cid: args.cid,
              conformer_type: args.conformer_type || '3d',
              properties: response.data,
            }, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get 3D conformers: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleAnalyzeStereochemistry(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid stereochemistry arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/property/AtomStereoCount,DefinedAtomStereoCount,BondStereoCount,DefinedBondStereoCount,IsomericSMILES/JSON`);

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({
              cid: args.cid,
              stereochemistry: response.data,
            }, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to analyze stereochemistry: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  private async handleGetCompoundProperties(args: any) {
    if (!isValidPropertiesArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid compound properties arguments');
    }

    try {
      const properties = args.properties || [
        'MolecularWeight', 'XLogP', 'TPSA', 'HBondDonorCount', 'HBondAcceptorCount',
        'RotatableBondCount', 'Complexity', 'HeavyAtomCount', 'Charge'
      ];

      const response = await this.apiClient.get(`/compound/cid/${args.cid}/property/${properties.join(',')}/JSON`);

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get compound properties: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  // Placeholder implementations for remaining methods (commented out - uncomment to implement)
  // private async handleCalculateDescriptors(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Descriptor calculation not yet implemented', args }, null, 2) }] };
  // }

  // private async handlePredictAdmetProperties(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'ADMET prediction not yet implemented', args }, null, 2) }] };
  // }

  // private async handleAssessDrugLikeness(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Drug-likeness assessment not yet implemented', args }, null, 2) }] };
  // }

  // private async handleAnalyzeMolecularComplexity(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Molecular complexity analysis not yet implemented', args }, null, 2) }] };
  // }

  // private async handleGetPharmacophoreFeatures(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Pharmacophore features not yet implemented', args }, null, 2) }] };
  // }

  // private async handleSearchBioassays(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Bioassay search not yet implemented', args }, null, 2) }] };
  // }

  private async handleGetAssayInfo(args: any) {
    try {
      const response = await this.apiClient.get(`/assay/aid/${args.aid}/JSON`);
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get assay info: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  // private async handleGetCompoundBioactivities(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Bioactivity search not yet implemented', args }, null, 2) }] };
  // }

  // private async handleSearchByTarget(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Target search not yet implemented', args }, null, 2) }] };
  // }

  // private async handleCompareActivityProfiles(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Activity profile comparison not yet implemented', args }, null, 2) }] };
  // }

  private async handleGetSafetyData(args: any) {
    if (!isValidCidArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid CID arguments');
    }

    try {
      const response = await this.apiClient.get(`/compound/cid/${args.cid}/classification/JSON`);
      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify(response.data, null, 2),
          },
        ],
      };
    } catch (error) {
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get safety data: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  // private async handleGetToxicityInfo(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Toxicity info not yet implemented', args }, null, 2) }] };
  // }

  // private async handleAssessEnvironmentalFate(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Environmental fate assessment not yet implemented', args }, null, 2) }] };
  // }

  // private async handleGetRegulatoryInfo(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Regulatory info not yet implemented', args }, null, 2) }] };
  // }

  // private async handleGetExternalReferences(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'External references not yet implemented', args }, null, 2) }] };
  // }

  // private async handleSearchPatents(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Patent search not yet implemented', args }, null, 2) }] };
  // }

  // private async handleGetLiteratureReferences(args: any) {
  //   return { content: [{ type: 'text', text: JSON.stringify({ message: 'Literature references not yet implemented', args }, null, 2) }] };
  // }

  private async handleBatchCompoundLookup(args: any) {
    if (!isValidBatchArgs(args)) {
      throw new McpError(ErrorCode.InvalidParams, 'Invalid batch arguments');
    }

    try {
      const results = [];
      for (const cid of args.cids.slice(0, 10)) {
        try {
          const response = await this.apiClient.get(`/compound/cid/${cid}/property/MolecularWeight,CanonicalSMILES,IUPACName/JSON`);
          results.push({ cid, data: response.data, success: true });
        } catch (error) {
          results.push({ cid, error: error instanceof Error ? error.message : 'Unknown error', success: false });
        }
      }

      return { content: [{ type: 'text', text: JSON.stringify({ batch_results: results }, null, 2) }] };
    } catch (error) {
      throw new McpError(ErrorCode.InternalError, `Batch lookup failed: ${error instanceof Error ? error.message : 'Unknown error'}`);
    }
  }

  private async handleGetPatentIds(args: any) {
    try {
      let cid: number | string;

      if (args.cid) {
        cid = args.cid;
      } else if (args.smiles) {
        // Resolve SMILES to CID first
        const cidResponse = await this.apiClient.get(`/compound/smiles/${encodeURIComponent(args.smiles)}/cids/JSON`);
        const cids = cidResponse.data?.IdentifierList?.CID;
        if (!cids || cids.length === 0) {
          return {
            content: [{ type: 'text', text: JSON.stringify({ message: 'No compound found for the given SMILES', smiles: args.smiles }, null, 2) }],
          };
        }
        cid = cids[0];
      } else {
        throw new McpError(ErrorCode.InvalidParams, 'Either "cid" or "smiles" is required for get_patent_ids');
      }

      const response = await this.apiClient.get(`/compound/cid/${cid}/xrefs/PatentID/JSON`);
      const patentIds = response.data?.InformationList?.Information?.[0]?.PatentID || [];

      return {
        content: [
          {
            type: 'text',
            text: JSON.stringify({
              cid,
              smiles: args.smiles || null,
              patent_count: patentIds.length,
              patent_ids: patentIds,
              patent_urls: patentIds.slice(0, 20).map((id: string) =>
                `https://patents.google.com/patent/${id}`
              ),
            }, null, 2),
          },
        ],
      };
    } catch (error) {
      if (error instanceof McpError) throw error;
      throw new McpError(
        ErrorCode.InternalError,
        `Failed to get patent IDs: ${error instanceof Error ? error.message : 'Unknown error'}`
      );
    }
  }

  async run() {
    const transport = new StdioServerTransport();
    await this.server.connect(transport);
    console.error('PubChem MCP server running on stdio');
  }
}

const server = new PubChemServer();
server.run().catch(console.error);

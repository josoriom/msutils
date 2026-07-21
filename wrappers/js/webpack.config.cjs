"use strict";
const path = require("node:path");
const fs = require("node:fs");
const webpack = require("webpack");

const version = require("./package.json").version;
const wasmPath = path.resolve(__dirname, "..", "..", "artifacts", version, "wasm", "quantion.wasm");

if (!fs.existsSync(wasmPath)) {
  throw new Error(`quantion: ${wasmPath} is missing. Run 'make wasm' at the repo root first.`);
}

const wasmDataUrl =
  "data:application/wasm;base64," + fs.readFileSync(wasmPath).toString("base64");

module.exports = {
  entry: "./src/index-wasm.ts",
  mode: "production",
  module: {
    rules: [
      {
        test: /\.ts$/,
        use: {
          loader: "ts-loader",
          options: {
            configFile: path.resolve(__dirname, "tsconfig.esm.json"),
            transpileOnly: true,
          },
        },
        exclude: /node_modules/,
      },
    ],
  },
  resolve: {
    extensions: [".ts", ".js"],
    fallback: {
      fs: false,
      path: false,
      url: false,
      util: false,
    },
  },
  output: {
    path: path.resolve(__dirname, "dist"),
    filename: "quantion.js",
    library: "quantion",
    libraryTarget: "umd",
    globalObject: "this",
  },
  plugins: [
    new webpack.DefinePlugin({
      __INLINE__: "true",
      __WASM_DATA_URL__: JSON.stringify(wasmDataUrl),
    }),
    new webpack.IgnorePlugin({ resourceRegExp: /^node:/ }),
  ],
  performance: { maxAssetSize: 2_000_000, maxEntrypointSize: 2_000_000 },
};

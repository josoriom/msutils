import { ContainmentCache } from "./containmentCache";

type FreeCell<H> = { handle: H | null; free: (h: H) => void };

export interface HandleEntry {
  cache: ContainmentCache;
}

type HandleMetadata = {
  cache: ContainmentCache;
  pending_lock: Promise<void>;
};

class HandleRegistry<H> {
  private readonly registry: FinalizationRegistry<FreeCell<H>> | null;
  private readonly metadata: Map<unknown, HandleMetadata>;

  constructor() {
    this.metadata = new Map();
    this.registry =
      typeof FinalizationRegistry !== "undefined"
        ? new FinalizationRegistry<FreeCell<H>>((cell) => {
            const h = cell.handle;
            if (h == null) return;
            cell.handle = null;
            try {
              cell.free(h);
            } catch (e: unknown) {
              console.warn("Error in HandleRegistry:", e);
            }
          })
        : null;
  }

  register(owner: object, handle: H, free: (h: H) => void): FreeCell<H> {
    const cell: FreeCell<H> = { handle, free };
    this.registry?.register(owner, cell, owner);
    this.metadata.set(handle, {
      cache: new ContainmentCache(),
      pending_lock: Promise.resolve(),
    });
    return cell;
  }

  release(owner: object, cell: FreeCell<H>): void {
    const h = cell.handle;
    if (h == null) return;
    cell.free(h);
    cell.handle = null;
    this.metadata.delete(h);
    this.registry?.unregister(owner);
  }

  get_entry(handle: H): HandleEntry {
    if (!this.metadata.has(handle)) {
      this.metadata.set(handle, {
        cache: new ContainmentCache(),
        pending_lock: Promise.resolve(),
      });
    }
    const meta = this.metadata.get(handle)!;
    return {
      cache: meta.cache,
    };
  }

  async with_lock<T>(
    handle: H,
    fn: (entry: HandleEntry) => Promise<T>,
  ): Promise<T> {
    const meta =
      this.metadata.get(handle) ||
      (() => {
        const m: HandleMetadata = {
          cache: new ContainmentCache(),
          pending_lock: Promise.resolve(),
        };
        this.metadata.set(handle, m);
        return m;
      })();

    const prev_lock = meta.pending_lock;

    let resolve_lock: () => void;
    const new_lock = new Promise<void>((resolve) => {
      resolve_lock = resolve;
    });
    meta.pending_lock = new_lock;

    await prev_lock;

    try {
      return await fn({
        cache: meta.cache,
      });
    } finally {
      resolve_lock!();
    }
  }
}

export const fileHandleRegistry = new HandleRegistry<unknown>();
